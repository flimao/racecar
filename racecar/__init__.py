#!/usr/env/python3
# -*- coding: utf-8 -*-

"""
racecar: This module is a framework with which to analyse and simulate racecar behavior.

By LmnICE
"""

__all__ = [ '__author__',
            'racecar',
            'series']

__author__ = 'Limão <dev@lmnice.me>'

from pint import UnitRegistry

u = UnitRegistry()

airdensity = 1.225 * u('kg/(m**3)')