#!/usr/env/python3
# -*- coding: utf-8 -*-

import csv
import pint

import racecar.series as Series

from . import ureg

class Racecar:
    """
    Class that describes a racecar, with its many systems (engine, suspension, etc),
    """

    def __init__(self, engine = None, trans = None, massd = None):

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
        pass


class Engine:
    """
    Class that describes the engine
    """

    def __init__(self, racecar, mass =0 * ureg.kilogram,
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

        if torque_data and torque_units:
            self.import_torque_data(csv_file = torque_data,
                                    units = torque_units)

    def torque(self, x = None):

        if x is None:
            return self.max_torque
        else:
            return self.torque_data(x)

    def import_torque_data(self, csv_file,
                           units = ('rpm', 'N.m')):
        # units specified in pint format

        global ureg

        self.torque_data = Series.Series(series = csv_file,
                                      units = units)

        self.max_torque = max(self.torque_data, key=lambda x: x[1])[1]

        max_hp = ureg('horsepower')
        for rpm, tq in self.torque_data:
            hp = (rpm * tq).to('horsepower')
            if hp > max_hp:
                max_hp = hp

        self.max_power = max_hp


class Transmission:
    """
    Class that describes the transmission
    """

    def __init__(self, racecar, ratios,
                 mass =0 * ureg.kilogram):

        # racecar is a Racecar object
        self.racecar = racecar
        racecar.trans = self

        # ratios is an iterable containing the final drive gearing ratio as first element,
        # then the gearing ratios for each gear

        self.__ratios = None
        self.ratios = ratios

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
