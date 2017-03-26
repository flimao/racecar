#!/usr/env/python3
# -*- coding: utf-8 -*-

import csv
import pint
import matplotlib.pyplot as plt
import numpy as np
import re

import racecar.series as Series

from . import ureg, airdensity

# define acceleration in G
ureg.define('gravity = 9.80665 m/(s**2) = G')

# define other non-Imperial, non-SI power units (from Wikipedia)
ureg.define('cavalo-vapor = 735 N.m/s = cv') # portuguese
ureg.define('Pferdestarke = 1 cv = PS') # german
ureg.define('cheval = 1 cv = ch') # french
ureg.define('paardenkracht = 1 cv = pk') # dutch
ureg.define('Лошадиная сила = 1 cv = лс') # russian
ureg.define('hästkraft = 1 cv = hk') # swedish
ureg.define('hevosvoima = 1 cv = hv') # finnish
ureg.define('lóerő = 1 cv = LE') # hungarian
ureg.define('cal-putere = 1 cv = CP') # romanian

class Racecar:
    """
    Class that describes a racecar, with its many systems (engine, suspension, etc),
    """

    def __init__(self, engine = None, trans = None, massd = None, tires = None,
                 body = None, mechloss = 0.15):

        global ureg

        # building a racecar with the following attributes
        if engine is None:
            self.engine = Engine(racecar = self)
        else:
            self.engine = engine

        if trans is None:
            self.trans = Transmission(racecar = self, ratios = [ 1, 1 ])
        else:
            self.trans = trans

        if massd is None:
            self.mass = MassDistribution(racecar = self,
                                         curb_mass = ureg('1300 kg'),
                                         length = ureg('4.644 m'),
                                         width = ureg('1.778 m'),
                                         height = ureg('1.473 m'))
        else:
            self.mass = massd

        if tires is None:
            self.tires = Tires(racecar = self,
                               spec = '225/45R17')

        if body is None:
            self.body = Body(racecar = self,
                             cx = 0.28, frontal_area = 1.5)

        self.mechloss = mechloss

    def wheel_torque(self, v_or_rpm, gear, mechloss = None):
        """
        Calculate total torque at the wheels.
        First parameter can be either speed or rpm. If it's speed, it's
            first converted into rpm given the gear.
        """

        if v_or_rpm.dimensionality == ureg('km/hr').dimensionality:
            rpm = self.rpm_from_speed_gear(v = v_or_rpm, gear = gear)
        else:
            rpm = v_or_rpm

        mechloss = mechloss or self.mechloss
        return self.engine.torque(rpm) * self.trans.ratio(gear) * (1-mechloss)

    def speed_from_rpm_gear(self, rpm, gear, unit = 'km/hr'):
        """
        Calculate speed given rpm and gear
        """
        speed = (rpm / self.trans.ratio(gear)) * self.tires.driven.fulld / 2

        return speed.to(ureg(unit))

    def acceleration_from_rpm_gear(self, rpm, gear,
                                   resistance = True, unit = 'G'):
        """
        Calculate how much acceleration can the engine provide at the given rpm
        and gear.

        If resistance is True, this subtract from wheel force the air and
        rolling resistances.
        """

        wheel_torque = self.wheel_torque(rpm, gear)
        wheel_radius = self.tires.driven.fulld/2
        wheel_force = wheel_torque / wheel_radius

        resist = 0
        if resistance:
            v = self.speed_from_rpm_gear(rpm = rpm, gear = gear)
            resist = self.body.resistance(v = v)
            resist += self.tires.rolling_resistance

        accel = (wheel_force - resist) / self.mass.mass

        return accel.to(ureg(unit))

    def weight_to_power(self):
        return self.mass.curb / self.engine.max_power

    def max_accel_at_gear(self, gear = 1, rpm_incr = 20 * ureg.rpm,
                          unit = 'G'):

        max_accel = 0 * ureg(unit).units
        max_accel_rpm = 0 * rpm_incr.units

        for rpm in np.linspace(self.engine.torque_data[0][0].magnitude,
                               self.engine.torque_data[-1][0].magnitude,
                               rpm_incr.magnitude) * rpm_incr.units:

            accel = self.acceleration_from_rpm_gear(rpm = rpm,
                                                           gear = gear)
            if accel > max_accel:
                max_accel = accel
                max_accel_rpm = rpm

        return [ max_accel_rpm, max_accel.to(unit) ]

    def min_rpm_tireslip(self, gear, rpm_incr = 20 * ureg.rpm, unit = 'G'):
        """
        Calculate the minimum RPM at selected gear that makes tires slip on launch
        """
        min_accel = self.tires.driven.max_accel

        for rpm in np.linspace(self.engine.torque_data[0][0].magnitude,
                               self.engine.torque_data[-1][0].magnitude,
                               rpm_incr.magnitude) * rpm_incr.units:

            accel = self.acceleration_from_rpm_gear(rpm=rpm,
                                                    gear=gear)
            if accel > min_accel:
                return [ rpm, accel.to(unit) ]

    def rpm_from_speed_gear(self, v, gear, unit = 'rpm'):
        wheel_speed = v / (self.tires.driven.fulld / 2)
        rpm = wheel_speed * self.trans.ratio(gear)

        return rpm.to(unit)


    def best_gear_at_speed(self, v):
        """
        Calculate the best gear to be in at a given speed. This will be the
        gear that transmits the most torque to the wheels
        """

        max_wtq = ureg('0 N.m')
        best_gear = 1

        for g in range(self.trans.ngears):

            rpm = self.rpm_from_speed_gear(v, gear = g + 1)
            if rpm > self.engine.redline:
                continue
            try:
                wheel_torque = self.wheel_torque(v_or_rpm = rpm, gear = g + 1)

            except ValueError: # RPM out of range
                wheel_torque = 0 * max_wtq.units

            if wheel_torque > max_wtq:
                best_gear = g + 1
                max_wtq = wheel_torque

        return best_gear

    def shiftpoints(self, v_incr = ureg('1 km/hr')):
        """
        Calculates the best shiftpoints for each gear
        Returns a list with each shiftpoint
        """

        max_theo_speed = self.speed_from_rpm_gear(rpm = self.engine.redline,
                                             gear = self.trans.ngears)

        shiftpoints = []

        n = max_theo_speed.magnitude / v_incr.magnitude

        for s in np.linspace(0,
                             max_theo_speed.magnitude,
                             n) * v_incr.units:
            best_gear_0 = self.best_gear_at_speed(v = s)
            best_gear_1 = self.best_gear_at_speed(v = s + v_incr)

            if best_gear_0 < best_gear_1:
                rpm = self.rpm_from_speed_gear(v = s, gear = best_gear_0)

                shiftpoints.append(rpm)

        return shiftpoints

    def top_speed(self, v_incr = ureg('1 km/hr')):
        """
        Estimates the top speed
        """

        max_theo_speed = self.speed_from_rpm_gear(rpm = self.engine.redline,
                                             gear = self.trans.ngears)

        top_speed = ureg('0 km/hr')

        resultant_prev = ureg('10**5 kgf')

        n = max_theo_speed.magnitude / v_incr.magnitude

        # will be starting at 20 km/hr, because there is no torque at standstill
        for s in np.linspace(20, max_theo_speed.magnitude,
                             n) * v_incr.units:
            g = self.best_gear_at_speed(v = s)

            wheel_torque = self.wheel_torque(v_or_rpm = s, gear = g)
            wheel_force = wheel_torque / (self.tires.driven.fulld / 2)
            resistance = self.body.resistance(v = s)
            resistance += self.tires.rolling_resistance

            resultant = wheel_force - resistance

            if resultant_prev > ureg('0 N') and resultant < ureg('0 N'):
                 return s
            else:
                resultant_prev = resultant

    def go(self, destination, shiftpoints,
                     v0 = ureg('0 km/hr'),
                     launchrpm = None,
                     trans_shift_time = ureg('0 s'),
                     time_incr = ureg('0.005 s'),
                     verbose = False):

        g0 = self.best_gear_at_speed(v=v0)
        rpm0 = self.rpm_from_speed_gear(v = v0, gear = g0)
        try:
            launchrpm = launchrpm or self.min_rpm_tireslip(gear = g0)[0]
        except TypeError: # no launchrpm, engine can`t slip tires at v0
            launchrpm = rpm0

        dist = 0.0 * ureg.meters
        speed = v0
        total_time = 0.0 * ureg.s
        total_time_launching = 0.0 * ureg.s

        variables = [ dist,
                      speed,
                      total_time]

        print('Initial speed is {0:3.0~f}'.format(speed))
        print('Initial gear is {0:n}'.format(g0))
        print('    End of simulation happens at ', end='')

        if destination.dimensionality == dist.dimensionality:
            stop_crit_index = 0
            print('total distance of {0:6.2~f}'.format(destination))

        elif destination.dimensionality == speed.dimensionality:
            stop_crit_index = 1
            print('final speed of {0:3.0~f}'.format(destination))

        else:
            stop_crit_index = 3
            print('total time of {0:4.0~f}'.format(destination))

        radius = self.tires.driven.fulld / 2

        # first gear
        # launching in launchrpm, riding the clutch so max tire accel is achieved
        rpm = self.rpm_from_speed_gear(v=speed, gear=g0)

        if rpm < launchrpm and g0 == 1:
            launching = True
            print('Launching from {0:5.0~f}'.format(launchrpm))

        else:
            launching = False

        # Go!
        g = g0
        i = 0

        accel_tire = self.tires.driven.max_accel

        stop_crit = variables[stop_crit_index]

        while stop_crit <= destination:
            i += 1

            if launching:
                total_time_launching += time_incr
                rpm = launchrpm

            wheel_torque = self.wheel_torque(rpm, gear=g)
            wheel_force = wheel_torque / radius
            resistance = self.body.resistance(v = speed)
            resistance += self.tires.rolling_resistance
            resultant = wheel_force - resistance

            accel_engine = resultant / self.mass.mass
            accel = min(accel_engine, accel_tire)
            if accel_engine > accel_tire:
                tire_slip = '(T)'
            else:
                tire_slip = ''

            dist += speed * time_incr + (accel / 2) * time_incr**2
            total_time += time_incr
            speed += accel * time_incr

            variables = [ dist, speed, total_time ]

            stop_crit = variables[stop_crit_index]

            if divmod(i, 100)[1] == 0:
                print('Still running... stop criterion is at {0:4.0~f}'.format(
                    stop_crit
                ))

                if verbose:
                    print('   Speed: {0:3.0~f}, Accel: {1:1.3~f} {2:s}'.format(speed,
                                                    accel.to('G'), tire_slip))
                    print('   RPM: {0:5.0~f}, Time: {1:7.3~f}'.format(rpm, total_time))


            rpm = self.rpm_from_speed_gear(v = speed, gear = g)
            if rpm > shiftpoints[g-1]:
                g += 1

                if verbose:
                    print('Shifted to gear {0:n} @ {1:4.0~f}'.format(g, speed))

                rpm = self.rpm_from_speed_gear(v = speed, gear = g)
                total_time += trans_shift_time


            if rpm > launchrpm and launching:
                launching = False

                if verbose:
                    print('Clutch fully engaged @ {0:4.0~f}'.format(speed))

        print('\nEnd of simulation!')
        print('Stats:')
        print('    Total time {0:7.3~f}'.format(total_time))
        print('    Total distance {0:5.1~f}'.format(dist))
        print('    Final speed {0:4.1~f}'.format(speed))

        if verbose:
            print('-----------')
            print('    Time launching {0:4.1~f}'.format(total_time_launching))
            print('    Total shifts {0:n}'.format(g-g0))

class Engine:
    """
    Class that describes the engine
    """

    def __init__(self, racecar, mass =0 * ureg.kilogram,
                 idle = 800 * ureg.rpm,
                 redline = None,
                 torque_data = None, torque_units = None):

        # racecar is a Racecar object
        self.racecar = racecar
        racecar.engine = self

        # mass is a quantity with mass units
        self.mass = mass

        self.torque_data = None
        self.power_data = None

        self.max_torque = None
        self.max_power = None

        self.idle = idle
        self.redline = redline

        if torque_data and torque_units:
            self.import_torque_data(csv_file = torque_data,
                                    units = torque_units)

    def torque(self, x = None):

        if x is None:
            return self.max_torque
        elif x not in self.torque_data:
            raise ValueError('{0:5.0~f} outside of torque data range.'.format(x))
        else:
            return self.torque_data(x)

    def power(self, x = None):
        if x is None:
            return self.max_power
        else:
            return (self.torque_data(x) * x).to('cv')

    def import_torque_data(self, csv_file,
                           units = ('rpm', 'N.m')):
        # units specified in pint format

        global ureg

        self.torque_data = Series.Series(series = csv_file,
                                      units = units)

        max_torque = max(self.torque_data, key=lambda x: x[1])
        self.max_torque = max_torque[1]
        self.max_torque_rpm = max_torque[0]

        self.redline = self.torque_data[-1][0]

        max_hp = ureg('cv')
        max_hp_rpm = self.idle
        for rpm, tq in self.torque_data:
            hp = (rpm * tq)
            if hp > max_hp:
                max_hp = hp
                max_hp_rpm = rpm

        self.max_power = max_hp.to('cv')
        self.max_power_rpm = max_hp_rpm

    def plot(self, power = True):
        rpms = []
        tqs = []
        powers = []

        for dp in self.torque_data:
            rpms.append(dp[0].magnitude)
            tqs.append(dp[1].to('kgf.m').magnitude)
            if power:
                pwr = dp[0]*dp[1]
                powers.append(pwr.to('cv').magnitude)

        ax = plt.subplot(111)
        ax2 = ax.twinx()
        ax.plot(rpms, tqs)
        ax.set_ylabel('Torque (m-kgf)')
        ax.set_xlabel('RPM')

        if power:
            ax2.plot(rpms, powers, 'r')
            ax2.set_ylabel('Power (cv)')

        plt.show()


class Transmission:
    """
    Class that describes the transmission
    """

    def __init__(self, racecar, ratios, driven = 'FWD',
                 mass =0 * ureg.kilogram):

        # racecar is a Racecar object
        self.racecar = racecar
        racecar.trans = self

        # ratios is an iterable containing the final drive gearing ratio as first element,
        # then the gearing ratios for each gear

        self.__ratios = None
        self.ratios = ratios

        # driven axle
        self.driven = driven

        # mass is a quantity with mass units
        self.mass = mass

    def ratio(self, gear):
        """"
        Returns the final gear ratio, including the diff ratio for the specified gear
        """
        try:
            return self.ratios[0] * self.ratios[gear]

        except IndexError:
            gear_mod = divmod(gear, 10)
            if gear_mod[0] == 1:
                gear_ord = 'th'
            else:
                if gear_mod[1] == 1: gear_ord = 'st'
                elif gear_mod[1] == 2: gear_ord = 'nd'
                elif gear_mod[1] == 3: gear_ord = 'rd'
                else: gear_ord = 'th'

            raise IndexError("There is no {0:n}{1} gear.".format(gear, gear_ord))

    @property
    def ratios(self):
        return self.__ratios

    @ratios.setter
    def ratios(self, ratio_list):

        try:
            ratio_list = [ float(x) for x in ratio_list ]

        except TypeError:
            raise

        if len(ratio_list) < 2: # must have at least the final diff and 1 gear
            raise IndexError("Must have at least the diff and 1 gear.")

        self.__ratios =  ratio_list
        self.ngears = len(self.__ratios) - 1

    def __getitem__(self, index):
        return self.ratios[index]

    def __len__(self):
        return self.ratios

    def __iter__(self):
        return iter(self.ratios)



class MassDistribution:
    """
    Class that describes how the racecar's mass is distributed along the frame
    """

    def __init__(self, racecar, curb_mass, length, width, height):

        # racecar is a Racecar object
        self.racecar = racecar
        racecar.mass = self

        # length, width, height is a quantity in length units
        self.length = length
        self.width = width
        self.height = height

        # mass is a quantity with mass units
        self.curb = curb_mass
        self.mass = self.curb + ureg('70 kg') + ureg('288.45 kg')

        self.frame_linear_density = (curb_mass - self.racecar.engine.mass -
                                    self.racecar.trans.mass) / self.length

class Tire:
    """
    Class that describes a tire (wheel diameter, tread width, etc)
    Accepts a string parameter (spec) which is converted to aforementioned attributes
    """

    def __init__(self, spec):

        self.load_index_range = {}
        self.speed_rating_range = {}

        # wheel diameter
        self.wheeld = None

        # construction (radial, diagonal, etc)
        self.construction = None

        # rubber height
        self.aspectratio = None
        self.rubberheight = None

        # rubber width
        self.width = None

        # application (passenger, light truck, etc)
        self.application = None

        # load rating
        self.load_index = None
        self.max_load = None

        # speed rating
        self.speed_rating = None
        self.max_speed = None

        # full diameter (wheel + rubber)
        self.fulld = None

        # max loads - units of acceleration
        self.max_accel = None
        self.max_brake = None
        self.max_lateral_load = None

        self.rolling_resistance = 0 # ureg('19 lbf')

        self.spec = spec


    @property
    def spec(self):
        return self.__spec

    @spec.setter
    def spec(self, spect):
        """
        Reads a tire spec and extracts the attributes
        This spec is modified to contain tire performance. After the regular
           tire spec, add the following:
           MA:N.NN    - max acceleration in Gs (optional, default = 1G)
           MB:N.NN    - max deceleration in Gs (optional, equal to MA if missing)
           ML:N.NN    - max lateral loads in Gs (optional, default = 1G)
        """

        global ureg


        default_MA = 1 * ureg.G
        default_ML = 1 * ureg.G

        re_str = r'(P|LT|ST|T)?(\d{2,3})(?:\/|-)(\d{2})(R|B|D|-)(\d{1,2})'
        re_str += r'(?:\s([0-9]{2,3})(\(?(?:A[1-8]|[B-Z])\)?))?'
        re_str += r'(?:\s(?:MA:(\d\.\d\d)))?'
        re_str += r'(?:\s(?:MB:(\d\.\d\d)))?'
        re_str += r'(?:\s(?:ML:(\d\.\d\d)))?'

        spec_re = re.match(re_str, spect)

        comparison = re.match(r'.*', 'comparison')

        if type(spec_re) == type(comparison):
            self.__spec = spect

            # application
            self.application = spec_re.group(1)

            # rubber width
            self.width = float(spec_re.group(2)) * ureg.millimeters

            # rubber height
            self.aspectratio = float(spec_re.group(3))/100
            self.rubberheight = self.aspectratio * self.width

            # construction
            self.construction = spec_re.group(4) or 'R'

            # wheel diameter
            self.wheeld = float(spec_re.group(5)) * ureg.inches

            # full diameter
            self.fulld = self.wheeld + 2 * self.rubberheight

            # load index
            try:
                self.load_index = float(spec_re.group(6))
            except TypeError: # no load index information
                self.load_index = None

            # speed rating
            self.speed_rating = spec_re.group(7)

            # max accels - max acceleration
            try:
                self.max_accel = float(spec_re.group(8)) * ureg.G
            except TypeError: # no max accel information
                self.max_accel = default_MA

            # max accels - max brake
            try:
                self.max_brake = float(spec_re.group(9)) * ureg.G
            except TypeError:  # no max brake information
                self.max_brake = self.max_accel

            # max accels - max lateral load
            try:
                self.max_lateral_load = float(spec_re.group(10)) * ureg.G
            except TypeError:  # no max lateral load information
                self.max_lateral_load = default_ML

        else:
            raise TypeError("Invalid tire specification {0}".format(spect))

    @property
    def load_index(self):
        return self.__load_index

    @load_index.setter
    def load_index(self, LI):
        self.__load_index = LI
        self.max_load = self.load_index_range.get(LI, None)

    @property
    def speed_rating(self):
        return self.__load_index

    @speed_rating.setter
    def speed_rating(self, SR):
        self.__speed_rating = SR
        self.max_speed = self.speed_rating_range.get(SR, None)

class Tires:
    """
    Class that describes the racecar's tires
    """

    def __init__(self, racecar, *args, **kwargs):

        # racecar is a Racecar object
        self.racecar = racecar
        racecar.tires = self

        if len(args) > 0:
            self.front = Tire(spec = args[0])
            self.rear = self.front

        elif 'spec' in kwargs:
            self.front = Tire(spec = kwargs['spec'])
            self.rear = self.front

        else:
            self.front = Tire(spec = kwargs.get('front_spec'))
            self.rear = Tire(spec = kwargs.get('front_spec'))

        # driven wheels
        if self.racecar.trans.driven == 'FWD':
            self.driven = self.front

        elif self.racecar.trans.driven == 'RWD':
            self.driven = self.rear

        else:
            self.driven = [ self.front, self.rear ]

        self.rolling_resistance = 2 * self.front.rolling_resistance + \
            2* self.rear.rolling_resistance

class Body:
    """"
    Class that describes the body aerodynamics.
    """

    def __init__(self, racecar, cx, frontal_area):

        # racecar is a Racecar object
        self.racecar = racecar
        racecar.body = self

        self.cx = cx
        self.frontal_area = frontal_area

    def resistance(self, v, power = False):
        """
        Calculate air resistance, in units of Force or Power
        """
        global airdensity

        resistance = (1/2) * airdensity * v**2 * self.cx * self.frontal_area
        resistance.ito('kgf')

        if power:
            resistance *= v
            resistance.ito(ureg.cv)

        return resistance

    def calibrate_area(self, v, power, unit = 'm**2'):
        """
        Calculate frontal area from speed and power consumed at said speed
        """

        global airdensity

        area = 2 * power / (airdensity * v**3 * self.cx)

        self.frontal_area = area.to(ureg(unit).units)
        return self.frontal_area
