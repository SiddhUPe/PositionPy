# 1. TIME
# For these calculations, it is convenient to use Julian dates.
import datetime
from decimal import Decimal, getcontext
import math

# set the precision for Decimal module
getcontext().prec = 50 # good precision, fast results

def julian_date(date=datetime.datetime.now()):
    time = date.timestamp() * 1000
    tzoffset = date.utcoffset().total_seconds() // 60 if date.utcoffset() else 0
    return Decimal((time / 86400000) - (tzoffset / 1440) + 2440587.5)

J1970 = Decimal('2440588')
dayMs = Decimal('86400000')
# this is the inverse of the Julian Date function

def fromJulian(j):
    return datetime.datetime.fromtimestamp(float((Decimal(j) + Decimal('0.5') - J1970) * (dayMs / Decimal('1000.0'))))
