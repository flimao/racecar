#!/usr/env/python3
# -*- coding: utf-8 -*-

import csv

import racecar as mRC

from . import ureg

class Series:
    """
    Class that provides facilities to work with data series, such as torque data
    """

    def __init__(self, series, units = None):
        # series must be a list with two columns, one X-axis and one Y-axis
        # X-axis must be in ascending order
        # units must be a two item list = [ x_axis_unit, y_axis_unit ]
        #     where both units are strings (pint units) or a pint object
        # if series list is already made up of pint objects, no units are necessary
        # in this case, if units are not empty, the series will be unit-converted

        self.__series = None
        self.__units = None

        if hasattr(series[0][0], 'units'):
            if units is None:
                units = [ x.units for x in series[0][:] ]

            else: # convert the series to the specified units
                series_conv = []
                for dp in series:
                    dp_conv = [ ]
                    for v, u in zip(dp, units):
                        dp_conv.append(v.to(u))

                    series_conv.append(dp_conv)

        self.units = units
        self.series = series

    def get_data_point(self, x):
        """
        Given the value on the x axis, find the corresponding y-axis value
        linear interpolation

        x must be a Quantity
        """
        for i in range(len(self.series[:])-1):
            dp_0 = self.series[i]
            dp_1 = self.series[i+1]

            if x < dp_0[0] or x > dp_1[0]:
                continue

            return dp_0[1] + (x - dp_0[0]) * (dp_1[1] - dp_0[1]) / (dp_1[0] - dp_0[0])

    def to(self, units):
        """
        Convert units in place.
        units is a list with two items, same as the unit() setter
        """

        global ureg # pint's unit registry

        comparison = ureg.meter

        newus = []
        for newu, oldu in zip(units, self.units):

            if newu is None:
                newu = oldu

            if type(newu) == type(comparison): # units are a pint object
                newus.append(newu)

            elif type(newu) == type(1*comparison):
                newus.append(newu.units)

            else:
                newus.append(ureg(newu).units)

        for data_point in self.series:
            for axis_v, axis_newu in zip(data_point, newus):
                axis_v.ito(axis_newu)



    @property
    def series(self):
        return self.__series

    @series.setter
    def series(self, series):

        # series can be a list with two columns, or a string to a csv file

        global ureg
        Q_ = ureg.Quantity

        units = self.units

        if type(series) == type('string'): # try and open the csv
            with open(series, mode = 'r') as csv_file:
                series_iterator = csv.reader(csv_file)
                series_iterator = list(series_iterator)

        else: # series must be a list
            series_iterator = series

        try:
            s = []
            for xy in series_iterator:
                l = []

                for axis_v, axis_u in zip(xy, units):
                    axis_v = Q_(float(axis_v), axis_u)
                    l.append(axis_v)

                s.append(l)
        except:
            raise

        self.__series = s

    @property
    def units(self):
        return self.__units

    @units.setter
    def units(self, units):

        # units is a tuple with two items, corresponding to each axis unit
        # units can be specified in string format, using pint's spec,
        # or by passing a pint unit object.

        global ureg # pint's unit registry

        comparison = ureg.meter

        try:
            us = []
            for u in units:

                if type(u) == type(comparison): # units are a pint object
                    us.append(u)
                else:
                    us.append(ureg(u).units)

        except:
            raise

        if self.series is not None:
            new_s = self.convert_units(us)
            self.__series = new_s

        self.__units = us

    def __call__(self, x):
        return self.get_data_point(x)

    def __getitem__(self, index):
        return self.series[index]

    def __len__(self):
        return len(self.series)

    def __iter__(self):
        return iter(self.series)

    def __contains__(self, item):
        x = item[0]
        lim_x_1 = self.series[0]
        lim_x_2 = self.series[-1]
        if lim_x_1 < x < lim_x_2:
            return True
        elif lim_x_2 < x < lim_x_1:
            return True
        else:
            return False