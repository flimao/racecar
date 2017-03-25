#!/usr/env/python3
# -*- coding: utf-8 -*-

"""
Run the racecar analyses available in the module
"""

import csv
import pint

import racecar as mRC

import racecar.racecar as RC

ureg = mRC.ureg

csvf_jetta_unichip = './data/torque-jetta-unichip.csv'
csvf_jetta_oem = './data/torque-jetta-oem.csv'
csvf_units = [ 'rpm', 'N.m' ]

jetta_gears = [ 4.06,
                3.46,
                2.05,
                1.3,
                0.9,
                0.704,
                0.588]

#torque = mRC.Series(series = csvf, units = ['revolutions/minute', 'N.m'])

jetta_unichip = RC.Racecar()

jetta_oem = RC.Racecar()

engine = RC.Engine(racecar = jetta_unichip,
                   torque_data = csvf_jetta_unichip, torque_units = csvf_units)

trans = RC.Transmission(racecar = jetta_unichip,
                        ratios = jetta_gears)

engine.torque_data.to([None, ureg('kgf.m')])

massd = RC.MassDistribution(racecar = jetta_unichip,
                            curb_mass =1650 * ureg('kg'),
                            length = 4.2 * ureg.meters)

rpm = engine.torque_data[0][0].units