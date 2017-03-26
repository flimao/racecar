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

csvf_piggyback = './data/torque-piggyback.csv'
csvf_oem = './data/torque-oem.csv'
csvf_units = [ 'rpm', 'N.m' ]

gears = [ 4.06,
            3.46,
            2.05,
            1.3,
            0.9,
            0.91 * 3.14/4.06,
            0.76 * 3.14/4.06]

trans_shift_time = ureg('100 milliseconds')

mechloss = 0.17

tire_spec_oem = r'225/45R17 91V MA:0.45 ML:0.60'
tire_spec_piggyback = r'225/45R17 91V MA:0.53 ML:0.75'

curb = ureg('1375 kg')

l = ureg('4.644 m')
w = ureg('1.778 m')
h = ureg('1.473 m')

cx = 0.28
frontal_area = ureg('3.3557 m**2')

#torque = mRC.Series(series = csvf, units = ['revolutions/minute', 'N.m'])

piggyback = RC.Racecar(mechloss = mechloss)
engine = RC.Engine(racecar = piggyback,
                   torque_data = csvf_piggyback, torque_units = csvf_units)
trans = RC.Transmission(racecar = piggyback,
                        ratios = gears)
tires = RC.Tires(racecar = piggyback,
                 spec = tire_spec_piggyback)
massd = RC.MassDistribution(racecar = piggyback,
                            curb_mass = curb,
                            length = l,
                            width = w,
                            height = h)
body = RC.Body(racecar = piggyback,
               cx = cx,
               frontal_area = frontal_area) # 29 cv at 120 km\h

oem = RC.Racecar(mechloss = mechloss)
engine_oem = RC.Engine(racecar = oem,
                   torque_data = csvf_oem, torque_units = csvf_units)
trans_oem = RC.Transmission(racecar = oem,
                        ratios = gears)
tires_oem = RC.Tires(racecar = oem,
                 spec = tire_spec_oem)
massd_oem = RC.MassDistribution(racecar = oem,
                            curb_mass = curb,
                            length = l,
                            width = w,
                            height = h)
body_oem = RC.Body(racecar = oem,
               cx = cx,
               frontal_area = frontal_area) # 29 cv at 120 km\h


engine.torque_data.to([ None, ureg('kgf.m').units ])

rpm = engine.torque_data[0][0].units

v60 = ureg('60 km/hr')
v80 = ureg('80 km/hr')
v100 = ureg('100 km/hr')
v120 = ureg('120 km/hr')
dqt = ureg('1/4 mi')
d400 = ureg('400 m')
d1000 = ureg('1000 m')

shifts = piggyback.shiftpoints()

shifts_oem = oem.shiftpoints()

#piggyback.go(v100, shiftpoints = shifts,
#             trans_shift_time = trans_shift_time,
#             verbose = True)

#piggyback.top_speed()

#oem.go(v100, shiftpoints = shifts_oem,
#             trans_shift_time = trans_shift_time,
#             verbose = True)

oem.top_speed()
