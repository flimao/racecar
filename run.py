#!/usr/env/python3
# -*- coding: utf-8 -*-

"""
Run the racecar analyses available in the module
"""

import numpy as np
import sympy.solvers as sv

import pint

import racecar as mRC

import racecar.racecar as RC

from racecar import u, airdensity

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

trans_shift_time = u('100 milliseconds')

mechloss = 0.17

tire_spec_oem = r'225/45R17 91V MA:0.45 ML:0.60'
tire_spec_piggyback = r'225/45R17 91V MA:0.53 ML:0.75'

curb = u('1375 kg')

l = u('4.644 m')
w = u('1.778 m')
h = u('1.473 m')
wheelbase = u('2.651 m')
wheelbase_rear = None

engine_mass = u('400 kg')
trans_mass = u('200 kg')

pointmasses = { 'engine': 'front', 'transmission': 'transaxle',
                'driver': 'standard'}

cx = 0.28
frontal_area = u('3.3557 m**2')

#torque = mRC.Series(series = csvf, units = ['revolutions/minute', 'N.m'])

piggyback = RC.Racecar(mechloss = mechloss)
engine = RC.Engine(racecar = piggyback,
                   torque_data = csvf_oem, torque_units = csvf_units,
                   mass = engine_mass)
trans = RC.Transmission(racecar = piggyback,
                        ratios = gears,
                        mass = trans_mass)
tires = RC.Tires(racecar = piggyback,
                 spec = tire_spec_piggyback)
massd = RC.MassDistribution(racecar = piggyback,
                            curb_mass = curb, dims = (l, w, h),
                            wheelbase = wheelbase, wheelbase_rear = wheelbase_rear,
                            pointmasses = pointmasses)
body = RC.Body(racecar = piggyback,
               cx = cx,
               frontal_area = frontal_area) # 29 cv at 120 km\h

oem = RC.Racecar(mechloss = mechloss)
engine_oem = RC.Engine(racecar = oem,
                       torque_data = csvf_oem, torque_units = csvf_units,
                       mass=engine_mass)
trans_oem = RC.Transmission(racecar = oem,
                            ratios = gears,
                            mass = trans_mass)
tires_oem = RC.Tires(racecar = oem,
                 spec = tire_spec_oem)
massd_oem = RC.MassDistribution(racecar = oem,
                            curb_mass = curb, dims = (l, w, h),
                            wheelbase = wheelbase, wheelbase_rear = wheelbase_rear,
                            pointmasses = pointmasses)
body_oem = RC.Body(racecar = oem,
               cx = cx,
               frontal_area = frontal_area) # 29 cv at 120 km\h


engine.torque_data.to([None, u('kgf.m').units])

rpm = engine.torque_data[0][0].units

v60 = u('60 km/hr')
v80 = u('80 km/hr')
v100 = u('100 km/hr')
v120 = u('120 km/hr')
dqt = u('1/4 mi')
d400 = u('400 m')
d1000 = u('1000 m')

#shifts = piggyback.shiftpoints()

#shifts_oem = oem.shiftpoints()

#piggyback.go(v100, shiftpoints = shifts,
#             trans_shift_time = trans_shift_time,
#             verbose = True)

#piggyback.top_speed()

#oem.go(v100, shiftpoints = shifts_oem,
#             trans_shift_time = trans_shift_time,
#             verbose = True)

#oem.top_speed()

cg = (0.65, 0.48, 0.5)
cg = massd.cg

k_dura = u('(100 kgf) / (1 cm)')
k_mole = u('(1 kgf) / (1 cm)')

P = (massd.mass * u('1 G')).to('kgf')

dN00 = (1-cg[0])*(1-cg[1])
dN01 = (1-cg[0])*(cg[1])
dN10 = (cg[0])*(1-cg[1])
dN11 = (cg[0])*(cg[1])

N00_est = P*dN00
N01_est = P*dN01
N10_est = P*dN10
N11_est = P*dN11

f = (N00_est + N01_est + N10_est + N11_est) / P

N00_en = N00_est / f
N01_en = N01_est / f
N10_en = N10_est / f
N11_en = N11_est / f

s_en = N00_en + N01_en + N10_en + N11_en

eqs, p0, p, ax = massd._montarsist(engine_torque=u('40 kgf.m'), cg=cg, k = k_dura)

n00 = p[0] * u('kgf')
n01 = p[1] * u('kgf')
n10 = p[2] * u('kgf')
n11 = (P - n00 - n01 - n10).to('kgf')
s = n00 + n01 + n10 + n11
ntx, nty, ntz = p[3:6]
nabs = np.sqrt(ntx**2 + nty**2 + ntz**2)
nx, ny, nz = np.array((-ntx, nty, ntz)) / nabs
theta_roll = np.arccos(nz / np.sqrt(nz**2 + ny**2)) * (180 / np.pi) * (ny/abs(ny))
theta_pitch = np.arccos(nz / np.sqrt(nz**2 + nx**2)) * (180 / np.pi) * (nx/abs(nx))
print('   P = {0:~5.0f}'.format(P))
print('      Calculated | Estimated')
print(' N00 = {0:~5.0f} | {1:~5.0f}'.format(n00, N00_en))
print(' N01 = {0:~5.0f} | {1:~5.0f}'.format(n01, N01_en))
print(' N10 = {0:~5.0f} | {1:~5.0f}'.format(n10, N10_en))
print(' N11 = {0:~5.0f} | {1:~5.0f}'.format(n11, N11_en))
print('____________________________')
print('Soma = {0:~5.0f} | {1:~5.0f}'.format(s, s_en))
print('Coords: (nx = {0:5.3f}, ny = {1:5.3f}, nz = {2:5.3f})'.format(nx, ny, nz))
print('        |n| = {0:5.3f}'.format(nabs))
print('       Accel = {0:6.3f} G'.format(ax))
print(' Angulo Roll = {0:6.3f} deg'.format(theta_roll))
print('Angulo Pitch = {0:6.3f} deg'.format(theta_pitch))
