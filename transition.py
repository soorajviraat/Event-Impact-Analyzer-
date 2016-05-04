from __future__ import division
import numpy as np
from math import *
from pandas import DataFrame, date_range
from datetime import date


def timeindex(date1):
    return float(date1.year + (date1.month-1)/12 + (date1.day-1)/365)


def stepfunction(dblstart, dblend, dbllaunch, dbltimenow, dblinterval=1/12):
    case1 = dbltimenow[(dbltimenow + dblinterval) <= dbllaunch]
    case2 = dbltimenow[dbltimenow >= dbllaunch]
    step_function_case1 = DataFrame.from_records(np.zeros([1, case1.shape[0]]), index=None, exclude=None, columns=case1)
    step_function_case2 = DataFrame.from_records(np.zeros([1, case2.shape[0]]), index=None, exclude=None, columns=case2)
    step_function_case1[:] = dblstart * dblinterval
    step_function_case2[:] = dblend * dblinterval
    step_function = step_function_case1.add(step_function_case2, fill_value=0)
    return step_function


def monthly_date_generator(start_date, end_date=""):
    period = (end_date.year - start_date.year)*12 + (end_date.month - start_date.month) + 1
    index = date_range(start_date, periods=period, freq='M')
    daterange = index.year + (index.month-1)/12 + (index.day-1)/365
    return daterange


def s_shapedcurve(dmax, dtimetomax, dbllag, vlaunch, vtime, dinterval):

    dbllaunch = timeindex(vlaunch)      # Convert time and launch values to double precision values
    dbltime = vtime

    if dtimetomax < 0:
        dtimetomax = 0

    if dtimetomax == dbllag or dtimetomax == 0:
        s_curve = stepfunction(0, dmax, dbllaunch + dtimetomax, dbltime, dinterval)
    else:
        dbllower = np.maximum(dbltime, dbllaunch)               # Evaluate the integral
        dblupper = np.maximum(dbltime + dinterval, dbllaunch)
        dblcon = 1 / (1 - dbllag / dtimetomax) * log((1 / 0.15 - 1) / ((1 / 0.98 - 1) ** (dbllag / dtimetomax)), e)
        dblslo = 1 / dtimetomax * (log((1 / 0.98 - 1), e) - dblcon)
        c_lower = dblcon + dblslo * (dbllower - dbllaunch)
        c_upper = dblcon + dblslo * (dblupper - dbllaunch)
        ec_lower = np.exp(c_lower)
        ec_upper = np.exp(c_upper)
        dbllowerlimit = dmax / dblslo * (c_lower - np.log1p(ec_lower))
        dblupperlimit = dmax / dblslo * (c_upper - np.log1p(ec_upper))
        s_curve = (dblupperlimit - dbllowerlimit) / dinterval
    return s_curve


def rapidcurve(dmax, dtimetomax, vlaunch, vtime, dinterval):

    dbllaunch = timeindex(vlaunch)  # Convert time and launch values to double precision values
    dbltime = vtime

    if dtimetomax <= 0:
            rapid_curve = stepfunction(0, dmax, dbllaunch, dbltime, dinterval)

    else:
        dbllower = np.maximum(dbltime, dbllaunch)  # Evaluate the integral
        dblupper = np.maximum(dbltime + dinterval, dbllaunch)
        dblslo = 3.50655789731998 / dtimetomax  # -3.50655789731998 comes from the Log(1-0.98) term in the TMax equation
        dbllowerlimit = dmax * (dbllower + 1 / dblslo * np.exp(-(dblslo * (dbllower - dbllaunch))))
        dblupperlimit = dmax * (dblupper + 1 / dblslo * np.exp(-(dblslo * (dblupper - dbllaunch))))
        rapid_curve = (dblupperlimit - dbllowerlimit) / dinterval
    return rapid_curve


def transition(dtype, dmax, dtimetomax, vtimeinlag, vlaunch, vtime, dinterval):
    # force curve type limits
    if dtype < 0:
        dtype = 0
    elif dtype > 10:
        dtype = 10

    # test for pure S-shaped curve, use lag parameter
    if vtimeinlag == 0:
        dtimelag = dtimetomax / 4
    else:
        dtimelag = vtimeinlag

    # calculate curve type weight
    spercentscurve = 1 - (dtype / 10)
    spercentrapidcurve = 1 - spercentscurve

    if dtimetomax < 0:
                    transitionfinal = 0
    else:
        transitionslow = spercentscurve * s_shapedcurve(dmax, dtimetomax, dtimelag, vlaunch, vtime, dinterval)
        transitionfast = spercentrapidcurve * rapidcurve(dmax, dtimetomax, vlaunch, vtime, dinterval)
        transitionfinal = transitionslow + transitionfast
    return transitionfinal

k = monthly_date_generator(date(2011, 12, 31), date(2020, 12, 31))
array1 = np.array(k)
print(transition(0, 200, 1, 0, date(2011, 12, 31), array1, 1/12))
