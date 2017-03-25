#!/usr/env/python3
# -*- coding: utf-8 -*-

"""
Run the racecar analyses available in the module
"""

import csv
import pint
import racecar as mRC

import racecar.racecar as RC

g_units = mRC.g_units

csvf_jetta_unichip = './data/torque-jetta-unichip.csv'
csvf_jetta_oem = './data/torque-jetta-oem.csv'
csvf_units = [ 'revolution/minute', 'N.m' ]

#torque = mRC.Series(series = csvf, units = ['revolutions/minute', 'N.m'])

jetta_unichip = RC.Racecar()

jetta_oem = RC.Racecar()

engine = RC.Engine(racecar = jetta_unichip,
                   torque_data = csvf_jetta_unichip, torque_units = csvf_units)

engine.torque_data.to([ None, g_units('kgf.m')])

massd = RC.MassDistribution(racecar = racecar,
                            curb_mass = 1650 * g_units('kg'),
                            length = 4.2 * g_units.meters)

rpm = engine.torque_data[0][0].units