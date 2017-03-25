#!/usr/env/python3
# -*- coding: utf-8 -*-

import csv
import pint
import re

import racecar.series as Series

from . import ureg

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

    def __init__(self, engine = None, trans = None, massd = None, tires = None):

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
                                         curb_mass =1600 * ureg('kilogram'),
                                         length =3.0 * ureg.meters)
        else:
            self.mass = massd

        if tires is None:
            self.tires = Tires(racecar = self,
                               spec = '225/45R17')

    def speed_from_rpm_gear(self, rpm, gear, unit = 'km/hr'):
        speed = (rpm / self.trans.ratio(gear)) * self.tires.driven.fulld / 2

        return speed.to(ureg(unit))

    def accelerateTo(self, destination, shiftpoints,
                     trans_shift_time = 0 * ureg.s,
                     time_incr = 0.5 * ureg.s):
        dist = 0.0
        speed = 0.0
        total_time = 0.0




class Engine:
    """
    Class that describes the engine
    """

    def __init__(self, racecar, mass =0 * ureg.kilogram, idle = 800 * ureg.rpm,
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

        if torque_data and torque_units:
            self.import_torque_data(csv_file = torque_data,
                                    units = torque_units)

    def torque(self, x = None):

        if x is None:
            return self.max_torque
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

        self.max_torque = max(self.torque_data, key=lambda x: x[1])[1]

        max_hp = ureg('cv')
        for rpm, tq in self.torque_data:
            hp = (rpm * tq)
            if hp > max_hp:
                max_hp = hp

        self.max_power = max_hp.to('cv')


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

    def __init__(self, racecar, curb_mass, length):

        # racecar is a Racecar object
        self.racecar = racecar
        racecar.mass = self

        # length is a quantity in length units
        self.length = length

        # mass is a quantity with mass units
        self.mass = curb_mass - self.racecar.engine.mass - self.racecar.trans.mass

        self.frame_linear_density = (curb_mass - self.racecar.engine.mass -
                                    self.racecar.trans.mass) / self.length

class Tire:
    """
    Class that dscribes a tire (wheel diameter, tread width, etc)
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
        re_str += r'(?:\s([0-9]{2,3})(\((?:A[1-8]|[B-Z])\)?))?'
        re_str += r'(?:\s(?:MA:(\d\.\d\d))?)?'
        re_str += r'(?:\s(?:MB:(\d\.\d\d))?)?'
        re_str += r'(?:\s(?:ML:(\d\.\d\d))?)?'

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

