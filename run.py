#!/usr/env/python3
# -*- coding: utf-8 -*-

"""
Run the racecar analyses available in the module
"""

import csv
import pint

import racecar as mRC

import racecar.racecar as RC

from racecar import ureg, airdensity

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

jetta_oem = RC.Racecar(mechloss = 0.85)

engine = RC.Engine(racecar = jetta_unichip,
                   torque_data = csvf_jetta_unichip, torque_units = csvf_units)

trans = RC.Transmission(racecar = jetta_unichip,
                        ratios = jetta_gears)

tires = RC.Tires(racecar = jetta_unichip,
                 spec = '225/45R17 91V MA:1.00 ML:0.75',
                 max_accel = 1 * ureg.G,
                 max_brake = 1 * ureg.G,
                 max_lateral_load = 0.75 * ureg.G)

massd = RC.MassDistribution(racecar = jetta_unichip,
                            curb_mass =1350 * ureg('kg'),
                            length = 4.2 * ureg.meters)

body = RC.Body(racecar = jetta_unichip,
               cx = 0.28,
               frontal_area = 3.3557 * ureg.m**2) # 29 cv at 120 km\h

engine.torque_data.to([ None, ureg('kgf.m').units ])

rpm = engine.torque_data[0][0].units