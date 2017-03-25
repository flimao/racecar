#!/usr/env/python3
# -*- coding: utf-8 -*-

import csv
import pint
import racecar as mRC
import racecar.series as Series

g_units = mRC.g_units

class Racecar:
    """
    Class that describes a racecar, with its many systems (engine, suspension, etc),
    """

    def __init__(self, engine = None, trans = None, massd = None):

        global g_units

        # building a racecar with the following attributes
        if engine is None:
            self.engine = Engine(self, mass = 0 * g_units('kilogram'))
        else:
            self.engine = engine

        if trans is None:
            self.trans = Transmission(self, mass = 0 * g_units('kilogram'))
        else:
            self.trans = trans

        if massd is None:
            self.mass = MassDistribution(self, curb_mass = 1600 * g_units('kilogram'),
                                         length = 3.0 * g_units.meters)
        else:
            self.mass = massd
        pass


class Engine:
    """
    Class that describes the engine
    """

    def __init__(self, racecar, mass = 0 * g_units.kilogram,
                 torque_data = None, torque_units = None):

        # racecar is a Racecar object
        self.racecar = racecar

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
                           units = ('revolutions/minute', 'N.m')):
        # units specified in pint format

        global g_units

        self.torque_data = Series.Series(series = csv_file,
                                      units = units)

        self.max_torque = max(self.torque_data, key=lambda x: x[1])[1]

        max_hp = g_units('horsepower')
        for rpm, tq in self.torque_data:
            hp = (rpm * tq).to('horsepower')
            if hp > max_hp:
                max_hp = hp

        self.max_power = max_hp


class Transmission:
    """
    Class that describes the transmission
    """

    def __init__(self, racecar, mass = 0 * g_units.kilogram):

        # racecar is a Racecar object
        self.racecar = racecar

        # mass is a quantity with mass units
        self.mass = mass


class MassDistribution:
    """
    Class that describes how the racecar's mass is distributed along the frame
    """

    def __init__(self, racecar, curb_mass, length):

        # racecar is a Racecar object
        self.racecar = racecar

        # length is a quantity in length units
        self.length = length

        # mass is a quantity with mass units
        self.mass = curb_mass - self.racecar.engine.mass - self.racecar.trans.mass

        self.frame_linear_density = (curb_mass - self.racecar.engine.mass -
                                    self.racecar.trans.mass) / self.length
