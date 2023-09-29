"""This module implements correctly rounded interval arithmetics according to Kahan-Novoa-Ratz
    according to Ratschek H., Rokne J. New computer methods for global optimization. – Halsted Press, 1988. book.
"""

import decimal as dec
import numpy as np

def set_precision(prec):
    """
    Sets the precision as the number of significant decimal digits in mantissa
    
    Parameters
    ----------
    prec : integer
        The number of decimal places (should be between 1 and some reasonable value)
    """
    prev = dec.getcontext().prec
    dec.getcontext().prec = prec
    return prev
    
def _set_rounding_mode(rounding_mode):  
    prev = dec.getcontext().rounding
    dec.getcontext().rounding = rounding_mode
    return prev
    
def _set_rounding_mode_default():
    return _set_rounding_mode(dec.ROUND_HALF_EVEN)

def _set_rounding_mode_ceil():
    return _set_rounding_mode(dec.ROUND_CEILING)

def _set_rounding_mode_floor():
    return _set_rounding_mode(dec.ROUND_FLOOR)


def _my_mul(a, b):
    if (a == dec.Decimal('0') and b == dec.Decimal('Infinity')) or (a == dec.Decimal('0') and b == dec.Decimal('-Infinity')) or (b == dec.Decimal('0') and a == dec.Decimal('Infinity')) or (b == dec.Decimal('0') and a == dec.Decimal('-Infinity')):
        return dec.Decimal('0')
    else:
        return a * b
    
def _my_div(a, b):
    if b.is_infinite():
        if a.is_infinite():
            raise ValueError('Infinity by infinity division')
        else:
            return dec.Decimal('0')
    return a / b

class Interval:
    """Class for storing interval values and perform interval operations"""
    
    def _convert_to_interval(other):
        if type(other) == dec.Decimal:
            return Interval(other, other)
        elif np.issubdtype(type(other), np.integer) or np.issubdtype(type(other), np.floating):
            v = dec.Decimal(other)
            return Interval(v, v)
        else:
            return other
        

    def __init__(self, a : [float, int, dec.Decimal, np.int64, np.float64], b : [float, int, dec.Decimal, np.int64, np.float64]):
        """
        Constructor

        Parameters
        ----------
        a : interval's left end
        b : interval's right end

        """
        if (not (np.issubdtype(type(a), np.integer) | np.issubdtype(type(a), np.floating)) and not type(a) == dec.Decimal)\
        or (not (np.issubdtype(type(b), np.integer) | np.issubdtype(type(b), np.floating)) and not type(b) == dec.Decimal):
             raise TypeError("Interval constructor's arguments must be Decimal or numerical")
        if a > b:
            raise ValueError("Left bound must be less or equal than the Right bound of interval")
        self.a = dec.Decimal(a)
        self.b = dec.Decimal(b)
        

    def __neg__(self):
        return Interval(-self.b, -self.a)

        
    def __add__(self, other):
        nother = Interval._convert_to_interval(other)
        _set_rounding_mode_floor()
        a = self.a + nother.a
        _set_rounding_mode_ceil()
        b = self.b + nother.b
        return Interval(a, b)
    
    def __sub__(self, other):
        nother = Interval._convert_to_interval(other)
        _set_rounding_mode_floor()
        a = self.a - nother.b
        _set_rounding_mode_ceil()
        b = self.b - nother.a
        return Interval(a, b)
    
    def __mul__(self, other):
        nother = Interval._convert_to_interval(other)
        _set_rounding_mode_floor()
        aa = self.a * nother.a
        ab = self.a * nother.b
        ba = self.b * nother.a
        bb = self.b * nother.b
        a = min(aa, ab, ba, bb)
        _set_rounding_mode_ceil()        
        aa = self.a * nother.a
        ab = self.a * nother.b
        ba = self.b * nother.a
        bb = self.b * nother.b
        b = max(aa, ab, ba, bb)
        return Interval(a, b)
    
    def mul(self, other): #multiplication in case of infinity bounds of intervals
        nother = Interval._convert_to_interval(other)
        _set_rounding_mode_floor()
        aa = _my_mul(self.a, nother.a)
        ab = _my_mul(self.a, nother.b)
        ba = _my_mul(self.b, nother.a)
        bb = _my_mul(self.b, nother.b)
        a = min(aa, ab, ba, bb)
        _set_rounding_mode_ceil()        
        aa = _my_mul(self.a, nother.a)
        ab = _my_mul(self.a, nother.b)
        ba = _my_mul(self.b, nother.a)
        bb = _my_mul(self.b, nother.b)
        b = max(aa, ab, ba, bb)
        return Interval(a, b)
        
    def __pow__(self, other):
        if (np.issubdtype(type(other), np.integer)) and other > 0:
            _set_rounding_mode_floor()
            if other % 2 == 0:
                if self.a <= 0 and self.b >= 0:
                    a = dec.Decimal('0')
                    _set_rounding_mode_ceil()
                    if -self.a < self.b:
                        b = self.b ** other
                    else:
                        b = self.a ** other
                elif self.a > 0:
                    a = self.a ** other
                    _set_rounding_mode_ceil()
                    b = self.b ** other
                elif self.b < 0:
                    a = self.b ** other
                    _set_rounding_mode_ceil()
                    b = self.a ** other
            else:
                a = self.a ** other
                _set_rounding_mode_ceil()
                b = self.b ** other                   
            return Interval(a, b)
        else:
            if (np.issubdtype(type(other), np.integer)) and other == 0:
                return Interval(1, 1)
            raise TypeError("Interval constructor's arguments must be a positive integer")
        
    
        
    def __truediv__(self, other):
        vzero = dec.Decimal('0')
        nother = Interval._convert_to_interval(other)
        if nother.a == nother.b == vzero:
            if self.a <= vzero <= self.b:
                return Interval(dec.Decimal('-Infinity'), dec.Decimal('Infinity'))
            else:
                return None
        elif nother.a > vzero or nother.b < vzero:
            _set_rounding_mode_floor()
            ra = dec.Decimal('1')/nother.b
            _set_rounding_mode_ceil()
            rb = dec.Decimal('1')/nother.a
            return self.__mul__(Interval(ra, rb))
        elif self.a <= vzero <= self.b:
            return Interval(dec.Decimal('-Infinity'), dec.Decimal('Infinity'))
        elif nother.a == vzero:
            if self.b < vzero:
                _set_rounding_mode_ceil()
                rb = _my_div(self.b, nother.b)
                return Interval(dec.Decimal('-Infinity'), rb)
            elif self.a > vzero:
                _set_rounding_mode_floor()
                ra = _my_div(self.a, nother.b)
                return Interval(ra, dec.Decimal('Infinity'))
        elif nother.b == vzero:
            if self.b < vzero:
                _set_rounding_mode_floor()
                ra = _my_div(self.b, nother.a)
                return Interval(ra, dec.Decimal('Infinity'))
            elif self.a > vzero:
                _set_rounding_mode_ceil()
                rb = _my_div(self.a, nother.a)
                return Interval(dec.Decimal('-Infinity'), rb)
        else: #nother.a < 0 < nother.b:
            if self.b < vzero:
                _set_rounding_mode_ceil()
                ra = _my_div(self.b, nother.b)
                _set_rounding_mode_floor()
                rb = _my_div(self.b, nother.a)
            elif self.a > vzero:
                _set_rounding_mode_ceil()
                ra = _my_div(self.a, nother.a)
                _set_rounding_mode_floor()
                rb = _my_div(self.a, nother.b)
            return Interval(dec.Decimal('-Infinity'), dec.Decimal('Infinity'))
        #[Interval(dec.Decimal('-Infinity'), ra), Interval(rb, dec.Decimal('Infinity'))]      
    
    def __radd__(self, other):
        nother = Interval._convert_to_interval(other)
        return nother.__add__(self)

    
    def __rsub__(self, other):
        nother = Interval._convert_to_interval(other)
        return nother.__sub__(self)

    
    def __rmul__(self, other):
        nother = Interval._convert_to_interval(other)
        return nother.__mul__(self)
    
    
    def __rtruediv__(self, other):
        nother = Interval._convert_to_interval(other)
        return nother.__truediv__(self)
    
    def __repr__(self):
        return "Interval[" + str(self.a) + ", " + str(self.b) + "]"
    
    
    
    def __eq__(self, other):
        if (other == None):
            return False
        nother = Interval._convert_to_interval(other)
        return (self.a == nother.a) & (self.b == nother.b)
    
    def dissect(self, trashold):
        trsh = dec.Decimal(trashold)
        if (trsh <= self.b):
            if (self.a <= trsh):
                left = Interval(self.a, trsh)
                right = Interval(trsh, self.b)
                return [left, right]
            else:
                return [None, self]
        else:
            return [self, None]

#     def mid(self):
#         """
#         Returns the best approximation of the interval central point
#         """
#         return dec.Decimal('0.5') * (self.a + self.b)

#     def mid_int(self):
#         """
#         Returns the minimal interval enclosing the interval central point
#         """
#         _set_rounding_mode_floor()
#         a = dec.Decimal('0.5') * (self.a + self.b)
#         _set_rounding_mode_ceil()
#         b = dec.Decimal('0.5') * (self.a + self.b)
#         return Interval(a, b)
    
    
    
    
# Some utility function for working with intervals

def convert_to_interval(val):
    """
    Returns an interval [ival, ival]
    
    Parameters
    ----------
    val : Decimal, int, float
    """
    if (type(val) == dec.Decimal) or np.issubdtype(type(val), np.integer) or np.issubdtype(type(val), np.floating):
        return Interval(val, val)
    else:
        raise TypeError("Arguments must be instances of Decimal")
    

def mid(ival):
    """
    Returns the best approximation of the interval central point
    
    Parameters
    ----------
    ival : interval
    """
    _set_rounding_mode_default()
    return dec.Decimal('0.5') * (ival.a + ival.b)

def mid_interval(ival):
    """
    Returns the narrowest interval, containing the interval central point
    
    Parameters
    ----------
    ival : interval
    """
    _set_rounding_mode_floor()
    a = dec.Decimal('0.5') * (ival.a + ival.b)
    _set_rounding_mode_ceil()
    b = dec.Decimal('0.5') * (ival.a + ival.b)
    return Interval(a, b)

    
def wid(ival):
    """
    Returns the closest outer approximation of the interval's width
    
    Parameters
    ----------
    ival : interval
    """
    _set_rounding_mode_ceil()
    return ival.b - ival.a 
    


def intersect_intervals(ival1, ival2):
    """
    Returns the intersection of intervals (None in there is no intersection)
    
    Parameters
    ----------
    ival1 : 1st interval
    ival2 : 2nd interval
    """
    
    if ival1.b < ival2.a or ival2.b < ival1.a:
        return None
    else:
        a = max(ival1.a, ival2.a)
        b = min(ival1.b, ival2.b)
        return Interval(a, b)

def isin_1d(ival1, ival2):
    return (ival2.a <= ival1.a) & (ival1.b <= ival2.b)
    
def convex_hull(ival1, ival2):
    if ival1 == None:
        return ival2
    if ival2 == None:
        return ival1
    return Interval(min(ival1.a, ival2.a), max(ival1.b, ival2.b))

def mignitude(ival):
    if ival.a > 0:
        return ival.a
    if ival.b < 0:
        return -ival.b
    return dec.Decimal(0)


def _h_(x1, area, p, meps):
    if mignitude(area - x1) < meps:
        return p * (convex_hull(area, x1) ** (p - 1))
    else:
        return ((x1 ** p) - (area ** p)) / (x1 - area)

def _slope_convex_pow_(X, area, p, meps):
    x1 = X.a
    x1 = Interval(x1, x1)
    x2 = X.b
    x2 = Interval(x2, x2)

    ival1 = _h_(x1, area, p, meps)

    ival2 = _h_(x2, area, p, meps)

    return convex_hull(ival1, ival2)

def pow_slope(X, area_x, p):
    if (not np.issubdtype(type(p), np.integer)) or (p < 0):
        raise ValueError("power p of X must be non-negative integer")
    if p == 0:
        return Interval(0, 0)
    if p == 1:
        return Interval(1, 1)
    
    meps = 1e-8
    
    area = Interval._convert_to_interval(area_x)
    
    if p % 2 == 0:
        return _slope_convex_pow_(X, area, p, meps)

    ans = None

    half_space = X.dissect(0)
    half_space_area = area.dissect(0)

#     if (half_space_area[0] != None):
#         if (half_space[0] != None):
#             ans = convex_hull(_slope_convex_pow_(half_space[0], half_space_area[0], p, meps), ans)
#         if (half_space[1] != None):
#             ans = convex_hull(_slope_convex_pow_(half_space[1], half_space_area[0], p, meps), ans)
#
#     if (half_space_area[1] != None):
#         if (half_space[1] != None):
#             ans = convex_hull(_slope_convex_pow_(half_space[1], half_space_area[1], p, meps), ans)
#         if (half_space[0] != None):
#             ans = convex_hull(_slope_convex_pow_(half_space[0], half_space_area[1], p, meps), ans)
    
    if (half_space[0] != None):
        ans = convex_hull(_slope_convex_pow_(half_space[0], area, p, meps), ans)
        if (half_space_area[1] != None):
            x0 = Interval(0, 0)
            ans = convex_hull(_h_(x0, half_space_area[1], p, meps) / 2, ans)

    if (half_space[1] != None):
        ans = convex_hull(_slope_convex_pow_(half_space[1], area, p, meps), ans)
        if (half_space_area[0] != None):
            x0 = Interval(0, 0)
            ans = convex_hull(_h_(x0, half_space_area[0], p, meps) / 2, ans)
    return ans

    
def pow_slope_wide(X, area_x, p):
    if (not np.issubdtype(type(p), np.integer)) or (p < 0):
        raise ValueError("power p of X must be non-negative integer")
    if p == 0:
        return Interval(0, 0)
    if p == 1:
        return Interval(1, 1)

    area = Interval._convert_to_interval(area_x)
    
    return intersect_intervals(pow_slope(X, area_x, p), p * (convex_hull(X, area) ** (p - 1)))

def contains_zero_1d(X):
    return (X.a <= 0) & (0 <= X.b)






class IntervalVector:
    __array_ufunc__ = None
    
    def _convert_to_interval_vector(other, shape=None):
        if np.issubdtype(type(other), np.integer):
            if shape != None:
                return IntervalVector(np.full(shape, dec.Decimal(other)), np.full(shape, dec.Decimal(other)))
            else:
                raise ValueError("Number has no shape")
        if type(other) == Interval:
            if shape != None:
                return IntervalVector(np.full(shape, other.a), np.full(shape, other.b))
            else:
                raise ValueError("Number has no shape")

        if type(other) == IntervalVector:
            return other
        elif (type(other) == np.ndarray):
            if (type(other.flat[0]) in [int, float, dec.Decimal, np.int64, np.float64]):
                return IntervalVector(other, other)
            elif (type(other.flat[0]) == Interval):
                return IntervalVector.from_vector(other)
        else:
            return other
    
    def __init__(self, a : np.array, b : np.array):
        """
        Constructor

        Parameters
        ----------
        a : interval's left end
        b : interval's right end

        """
        if (a.shape != b.shape):
            raise ValueError(f"Bounds must have same shape {a.shape} vs {b.shape}")

        if np.any(a > b):
            raise ValueError("Left bound must be less or equal than the Right bound of interval")
        
        self.a = np.vectorize(lambda x: dec.Decimal(x))(a.astype(object))
        self.b = np.vectorize(lambda x: dec.Decimal(x))(b.astype(object))
        self.shape = self.a.shape
    
    def __len__(self):
        return len(self.a)
    
    def is_eq(self, other):
        return (self.a == other.a).all() & (self.b == other.b).all()
    
    def from_vector(other: np.ndarray):
        """
        Construct IntervalVector from usual array of Intervals
        """
        if (type(other) == IntervalVector):
            return other
        dist = Interval
        if np.vectorize(lambda x: type(x) in [dist, float, int])(other).all():
            return IntervalVector(np.vectorize(lambda x: Interval._convert_to_interval(x).a)(other), np.vectorize(lambda x: Interval._convert_to_interval(x).b)(other))
        else:
            raise TypeError('All elements must be Intervals')
    
    def to_vec(self):
        """
        Construct usual np.array of Intervals from IntervalVector
        """
        return np.apply_along_axis(lambda x: Interval(x[0], x[1]), 0, np.stack([self.a, self.b]))

    def __neg__(self):
        return IntervalVector(-self.b, -self.a)

        
    def __add__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        _set_rounding_mode_floor()
        a = self.a + nother.a
        _set_rounding_mode_ceil()
        b = self.b + nother.b
        return IntervalVector(a, b)
    
    def __sub__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        _set_rounding_mode_floor()
        a = self.a - nother.b
        _set_rounding_mode_ceil()
        b = self.b - nother.a
        return IntervalVector(a, b)
    
    def __mul__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        _set_rounding_mode_floor()
        aa = self.a * nother.a
        ab = self.a * nother.b
        ba = self.b * nother.a
        bb = self.b * nother.b
        a = np.stack([aa, ab, ba, bb]).min(axis=0)
        _set_rounding_mode_ceil()        
        aa = self.a * nother.a
        ab = self.a * nother.b
        ba = self.b * nother.a
        bb = self.b * nother.b
        b = np.stack([aa, ab, ba, bb]).max(axis=0)
        return IntervalVector(a, b)
    
    def sum(self, ax=None):
        _set_rounding_mode_floor()
        a = self.a.sum(axis=ax)
        _set_rounding_mode_ceil()
        b = self.b.sum(axis=ax)
        return Interval(a, b)

    
    def __matmul__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other)
        if (len(nother.shape) > 2) or (len(self.shape) > 2):
            raise ValueError("Matmul is only for matrices or 1d-arrays")
        
        sq0 = False
        sq1 = False
        if (len(nother.shape) == 1):
            nother = IntervalVector(nother.a[:, np.newaxis], nother.b[:, np.newaxis])
            sq1 = True
        if (len(self.shape) == 1):
            this = IntervalVector(self.a[np.newaxis, :], self.b[np.newaxis, :])
            sq0 = True
        else:
            this = self
        
        x_a = np.repeat(this.a, nother.shape[1], axis=0)
        x_b = np.repeat(this.b, nother.shape[1], axis=0)
        y_a = np.tile(nother.a.T, (this.shape[0], 1))
        y_b = np.tile(nother.b.T, (this.shape[0], 1))
        
        _set_rounding_mode_floor()
        aa = x_a * y_a
        ab = x_a * y_b
        ba = x_b * y_a
        bb = x_b * y_b
        a = np.stack([aa, ab, ba, bb]).min(axis=0)
        a = a.sum(axis=1).reshape(this.shape[0], nother.shape[1])

        _set_rounding_mode_ceil()        
        aa = x_a * y_a
        ab = x_a * y_b
        ba = x_b * y_a
        bb = x_b * y_b
        b = np.stack([aa, ab, ba, bb]).max(axis=0)
        b = b.sum(axis=1).reshape(this.shape[0], nother.shape[1])
        
        if sq0 and sq1:
            a = a.squeeze()
            b = b.squeeze()
            return Interval(a.item(), b.item())
        elif sq1:
            a = a.squeeze(axis=1)
            b = b.squeeze(axis=1)
        elif sq0:
            a = a.squeeze(axis=0)
            b = b.squeeze(axis=0)
        
        return IntervalVector(a, b)

    def __and__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        a = np.stack([self.a, nother.a]).max(axis=0)
        b = np.stack([self.b, nother.b]).min(axis=0)
        if (a <= b).all():
            return IntervalVector(a, b)
        else:
            return None
        
    def __pow__(self, other):
        return self.to_vec() ** other
    
    def __getitem__(self, i):
        if type(i) not in [slice, tuple]:
            a = self.a[i]
            b = self.b[i]
            if type(a) == dec.Decimal:
                return Interval(a, b)
            return IntervalVector(a, b)
        else:
            a = self.a[i]
            b = self.b[i]
            if len(a) == 0:
                return None
            return IntervalVector(self.a[i], self.b[i])
    
    def squeeze(self):
        return IntervalVector(self.a.squeeze(), self.b.squeeze())
    
    def item(self, *ind):
        return Interval(self.a[ind], self.b[ind])
    
    def M(self):
        return np.stack([self.a, -self.a, self.b, -self.b]).max()
    
    def convex_hull(self, other):
        if other == None:
            return self
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        a = np.stack([self.a, nother.a]).min(axis=0)
        b = np.stack([self.b, nother.b]).max(axis=0)
        return IntervalVector(a, b)

    def concat(self, other):
        if (other is None):
            return self
        if (len(self.shape) != 1) or (len(other.shape) != 1):
            raise ValueError("both concatenating objects must be 1d vectors")
        nother = IntervalVector._convert_to_interval_vector(other)
        a = np.concatenate([self.a, nother.a])
        b = np.concatenate([self.b, nother.b])
        return IntervalVector(a, b)
        
        
#     def __truediv__(self, other):
#         vzero = dec.Decimal('0')
#         nother = IntervalVector._convert_to_interval_vector(other, self.a.shape)
#         if nother.a == nother.b == vzero:
#             if self.a <= vzero <= self.b:
#                 return Interval(dec.Decimal('-Infinity'), dec.Decimal('Infinity'))
#             else:
#                 return None
#         elif nother.a > vzero or nother.b < vzero:
#             _set_rounding_mode_floor()
#             ra = dec.Decimal('1')/nother.b
#             _set_rounding_mode_ceil()
#             rb = dec.Decimal('1')/nother.a
#             return self.__mul__(Interval(ra, rb))
#         elif self.a <= vzero <= self.b:
#             return Interval(dec.Decimal('-Infinity'), dec.Decimal('Infinity'))
#         elif nother.a == vzero:
#             if self.b < vzero:
#                 _set_rounding_mode_ceil()
#                 rb = _my_div(self.b, nother.b)
#                 return Interval(dec.Decimal('-Infinity'), rb)
#             elif self.a > vzero:
#                 _set_rounding_mode_floor()
#                 ra = _my_div(self.a, nother.b)
#                 return Interval(ra, dec.Decimal('Infinity'))
#         elif nother.b == vzero:
#             if self.b < vzero:
#                 _set_rounding_mode_floor()
#                 ra = _my_div(self.b, nother.a)
#                 return Interval(ra, dec.Decimal('Infinity'))
#             elif self.a > vzero:
#                 _set_rounding_mode_ceil()
#                 rb = _my_div(self.a, nother.a)
#                 return Interval(dec.Decimal('-Infinity'), rb)
#         else: #nother.a < 0 < nother.b:
#             if self.b < vzero:
#                 _set_rounding_mode_ceil()
#                 ra = _my_div(self.b, nother.b)
#                 _set_rounding_mode_floor()
#                 rb = _my_div(self.b, nother.a)
#             elif self.a > vzero:
#                 _set_rounding_mode_ceil()
#                 ra = _my_div(self.a, nother.a)
#                 _set_rounding_mode_floor()
#                 rb = _my_div(self.a, nother.b)
#             return [Interval(dec.Decimal('-Infinity'), ra), Interval(rb, dec.Decimal('Infinity'))]
            
            
            
       
        
    
    def __radd__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        return nother.__add__(self)

    
    def __rsub__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        return nother.__sub__(self)

    
    def __rmul__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        return nother.__mul__(self)
    
    def __rmatmul__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        return nother.__matmul__(self)
    
    def __rand__(self, other):
        nother = IntervalVector._convert_to_interval_vector(other, self.shape)
        return nother.__and__(self)

    
    def dot_mid(self, rounding_mod=None):
        if (rounding_mod == "ceil"):
            _set_rounding_mode_ceil()
        elif (rounding_mod == "floor"):
            _set_rounding_mode_floor()
        elif (rounding_mod != None):
            raise ValueError("Parameter rounding_mod must be 'ceil', 'floor' or None")
        
        mid = (self.a + self.b) / 2

        if isin(mid, self):
            return np.vectorize(lambda x: float(x))(mid)
        else:
            return np.vectorize(lambda x: float(x))(self.a)
        
    def fast_mid(self): # fast but doesn't always guarantee, that mid is in interval if interwal.w -> machine epsilon
         return ((self.a + self.b) / 2).astype(float)
        
    
    def mid(self):
        _set_rounding_mode_floor()
        a = (self.a + self.b) / 2
        
        _set_rounding_mode_ceil()
        b = (self.a + self.b) / 2
        return IntervalVector(np.stack([a, self.a]).max(axis=0), np.stack([b, self.b]).min(axis=0))

    
#     def __rtruediv__(self, other):
#         nother = IntervalVector._convert_to_interval_vector(other, self.a.shape)
#         return nother.__truediv__(self)
    
    def __repr__(self):
        return str(self.to_vec())
    
    def w(self):
        _set_rounding_mode_ceil()
        return (self.b - self.a).max()

    def coordinate_w(self):
        _set_rounding_mode_ceil()
        return (self.b - self.a)
    
    def dissect(self, ind, trashold):
        trsh = dec.Decimal(trashold)
        if len(self.shape) != 1:
            raise ValueError("Dissected IntervalVector must be 1d-vector")
        if (trsh <= self.b[ind]):
            if (self.a[ind] <= trsh):
                left = IntervalVector(self.a, np.concatenate([np.append(self.b[:ind], trsh), self.b[ind + 1:]]))
                right = IntervalVector(np.concatenate([np.append(self.a[:ind], trsh), self.a[ind + 1:]]), self.b)
                return [left, right]
            else:
                return [None, self]
        else:
            return [self, None]


    
def interval_scalar_prod(one, other):
        y = IntervalVector._convert_to_interval_vector(other)
        x = IntervalVector._convert_to_interval_vector(one)
        return (x * y).sum()

def isin(one, other):
    this = IntervalVector._convert_to_interval_vector(one)
    nother = IntervalVector._convert_to_interval_vector(other)
    if (nother.a <= this.a).all() and (this.b <= nother.b).all():
        return True
    return False

def contains_zero(other):
    nother = IntervalVector._convert_to_interval_vector(other)
    return (nother.a <= 0).all() & (0 <= nother.b).all()

def Interval2IntervalVector(one, shape):
    return IntervalVector(np.full(shape, one.a), np.full(shape, one.b))

def IntervalExtension(f):
    return lambda x: IntervalVector.from_vector(f(x))


def HansenSlope(X_wide, area_x, slope, n):
#     X = X_wide.convex_hull(area_x)
    ans = []
    for i in range(n):
        X_new = X_wide[:i + 1].concat(area_x[i + 1:])
        ans += [slope(X_new, area_x, i)]
    return np.array(ans)

def is_max_rank_floor(Y, tol=1e-7):
    Y_ival = np.vectorize(lambda x: convert_to_interval(x))(Y)
    if (len(Y_ival.shape) != 2):
        raise ValueError("To check the rank matrix must be 2-dimentional")
    if (Y_ival.shape[0] < Y_ival.shape[1]):
        Y_ival = Y_ival.T
    n = Y_ival.shape[0]
    m = Y_ival.shape[1]
    for j in range(m):
        i = 0

        while (i < n) and contains_zero_1d(Y_ival[i][j]):
            i += 1
        if i == n:
            return Y, False
        col = Y_ival[:, j]
        for k in range(j + 1, m):
            Y_ival[:, k] = Y_ival[:, k] - col * (Y_ival[i][k] / Y_ival[i][j])
            Y_ival[i][k] = Interval(0, 0)
    Y_ival_saved = np.vectorize(lambda x: convert_to_interval(x))(Y)
    return Y_ival_saved, True




from IPython.display import display, Math

class Monomial:
    def __init__(self, n, p={}):
        for i in p:
            if (not np.issubdtype(type(p[i]), np.integer)) or (p[i] <= 0):
                raise ValueError("All coefficients of Monomial must be positive integers")
            if (not np.issubdtype(type(i), np.integer)) or (i < 0) or (i >= n):
                raise ValueError("All indecies of Monomial must be positive integers between 0 and n")
        self.p = p
        self.n = n
    
    def __call__(self, x):
        if len(x) != self.n:
            raise ValueError("Inittial number of arguments must be same with number of arguments in call function")
        
        ans = 1
        for i in self.p:
            ans = ans * (x[i] ** self.p[i])
        return ans
    
    def grad(self, x):
        if len(x) != self.n:
            raise ValueError("Inittial number of arguments must be same with number of arguments in call function")
        
        ans = np.zeros(self.n, dtype=object)
        
        left_cummul = {}
        right_cummul = {}
        ind = []

        last_res = 1
        for i in self.p:
            ind += [i]
            left_cummul[i] = last_res
            last_res = last_res * (x[i] ** self.p[i])
        
        last_res = 1
        for i in ind[::-1]:
            right_cummul[i] = last_res
            last_res = last_res * (x[i] ** self.p[i])

        for i in self.p:
            power = self.p[i]
            if power == 1:
                arg = 1
            else:
                arg = power * (x[i] ** (power - 1))
            ans[i] = arg * left_cummul[i] * right_cummul[i]
        return ans
    
    def exclude_call(self, x, j):
        if len(x) != self.n:
            raise ValueError("Inittial number of arguments must be same with number of arguments in call function")
        if (not np.issubdtype(type(j), np.integer)) or (j < 0) or (j >= self.n):
                raise ValueError("All indecies of Monomial must be positive integers between 0 and n")
        ans = 1
        for i in self.p:
            if i != j:
                ans = ans * (x[i] ** self.p[i])
        return ans
    
    def to_str(self):
        ans = ''
        for ind in self.p:
            ans += f'x_{ind}^{self.p[ind]}'
        return ans

class Polynomial:
    def __init__(self, n, l, coeff):
        if len(coeff) != len(l):
            raise ValueError("Number of coefficients must be equal to the Number of monomials")
        mono_list = []
        for m in l:
            mono_list += [Monomial(n, m)]
        self.n = n
        self.coeff = coeff
        self.l = mono_list
        self.k = len(coeff)
        
        self.map = [[] for i in range(n)]
        for j, mono in enumerate(self.l):
            for i in mono.p:
                self.map[i] += [j]

        
    def __call__(self, x):
        if len(x) != self.n:
            raise ValueError("Inittial number of arguments must be same with number of arguments in call function")
        ans = 0
        for i in range(self.k):
            ans += self.coeff[i] * self.l[i](x)
        return ans
    
    def grad(self, x):
        if len(x) != self.n:
            raise ValueError("Inittial number of arguments must be same with number of arguments in call function")
        
        ans = np.zeros(self.n, dtype=object)
        for i in range(self.k):
            ans += self.coeff[i] * self.l[i].grad(x)
        
        return ans

    def slope_1d_from_nd(self, X, dot_x, i): #only if dot_x in X
        ans = 0
        xi = X[i]
        dot_xi = dot_x[i]
        for mono in self.map[i]:
            ans += self.coeff[mono] * pow_slope(xi, dot_xi, self.l[mono].p[i]) * self.l[mono].exclude_call(X, i)
        return ans
    
    def slope_nd(self, X, dot_x):
        if len(X.shape) != 1 or len(X) != self.n:
            raise ValueError("input must be n-dimensional vector")
        if len(dot_x.shape) != 1 or len(dot_x) != self.n:
            raise ValueError("input must be n-dimensional vector")
        
        return HansenSlope(X, dot_x, self.slope_1d_from_nd, self.n)
        
    
    def to_str(self):
        ans = ''
        for i in range(self.k):
            if i == 0:
                ans += str(self.coeff[i]) + self.l[i].to_str()
            else:
                if self.coeff[i] < 0:
                    ans += '-' + str(-self.coeff[i]) + self.l[i].to_str()
                else:
                    ans += '+' + str(self.coeff[i]) + self.l[i].to_str()
        return ans
    
class Multinomial:
    def __init__(self, n, l):
        self.ax = []
        for mono, coeff in l:
            self.ax += [Polynomial(n, mono, coeff)]
        self.n = n
        self.m = len(l)
    
    def __call__(self, x):
        if len(x) != self.n:
            raise ValueError("Inittial number of arguments must be same with number of arguments in call function")
        ans = []
        for poly in self.ax:
            ans += [poly(x)]
        return np.array(ans, dtype=object)
    
    def grad(self, x):
        if len(x) != self.n:
            raise ValueError("Inittial number of arguments must be same with number of arguments in call function")
        ans = []
        for poly in self.ax:
            ans += [poly.grad(x)]
        return np.stack(ans)
    
    def slope(self, X, dot_x):
        dot_x = IntervalVector.from_vector(dot_x)
        ans = []
        for poly in self.ax:
            ans += [poly.slope_nd(X, dot_x)]
        ans = np.stack(ans)

        return IntervalVector.from_vector(ans)
    
    def print_self(self):
        print(f'argument dimentionality: {self.n}, multinomial dimentionality: {self.m}')
        for i, p in enumerate(self.ax):
            display(Math(f'y_{i} = ' + p.to_str()))
        print('\n')
        
        
# n - input dimention
# variables_num - number of variables in monomial
# max_degree - max power of one variable in x
def random_monomial(n, variables_num, max_degree=10, seed=None): # produce monomial with k powers of different arguments
    if seed != None:
        np.random.seed(seed)

    x = np.random.permutation(n)[:variables_num]
    y = np.random.randint(1, max_degree, size=variables_num).astype(object)
    d = {}
    for i, power in zip(x, y):
        d[i] = power
    return d

# n - input dimention
# monomials_num - number of non-constant monomials
# const - constant monomial itself (just number)
# max_degree - Monomial parameter - max power of one variable in x
# max_monomial_variables_num - max number of variables in every monomial
# max_coeff - max module of coefficients before Monomials 

def random_polynomial(n, monomials_num, const=3, max_degree=10, max_coeff=15, max_monomial_variables_num=5, seed=None): # produce k monomials + 1 const monomial
    if seed != None:
        np.random.seed(seed)

    coeff = np.random.randint(-max_coeff, max_coeff, size=monomials_num).astype(object)
    coeff = np.append(coeff, const)
    
    monomial_k = np.random.randint(1, max_monomial_variables_num, size=monomials_num).astype(object)
    monoms = []
    
    for i in range(monomials_num):
        monoms += [random_monomial(n, monomial_k[i], max_degree)]
    monoms += [{}]
    return monoms, coeff

# n - input dimention
# m - output dimention
# max_monomials_num - max number of non-constant monomials per Polynomial
# max_module_constant - max module of constant in Polynomials
# max_degree - Monomial parameter: max power of one variable in x
# max_monomial_variables_num - max number of variables in every monomial
# max_coeff - max module of coefficients before Monomials

#produces m random Polynomials
def random_multinomial(n, m, max_monomials_num=10, max_module_constant=15, max_degree=10, max_coeff=15, max_monomial_variables_num=5, seed=None):
    if seed != None:
        np.random.seed(seed)
    ans = []
    monomials_num = np.random.randint(0, max_monomials_num, size=m).astype(object)
    polynomial_constants = np.random.randint(-max_module_constant, max_module_constant, size=m).astype(object)
    
    for i in range(m):
        ans += [random_polynomial(n, monomials_num[i], const=polynomial_constants[i], max_degree=max_degree, max_coeff=max_coeff, max_monomial_variables_num=max_monomial_variables_num)]
    
    return ans
    
        
from tqdm import tqdm

def test_set(seed=1234, number=15):
    np.random.seed(seed)
    n = np.random.randint(1, 10, size=number)
    test = []
    for i in range(number):
        test += [Multinomial(n[i], random_multinomial(n[i], n[i], max_monomials_num=7, max_module_constant=7, max_degree=5, max_coeff=15, max_monomial_variables_num=5))]
    return test, n

def test_set_nonsquare(seed=1234, number=15):
    np.random.seed(seed)
    n = np.random.randint(1, 10, size=number)
    m = np.random.randint(1, 10, size=number)
    tmp = np.stack([n, m])
    n = tmp.max(axis=0)
    m = tmp.min(axis=0)
    test = []
    for i in range(number):
        test += [Multinomial(n[i], random_multinomial(n[i], m[i], max_monomials_num=7, max_module_constant=7, max_degree=5, max_coeff=15, max_monomial_variables_num=5))]
    return test, n, m

def _test_(algo, number=15, **params):
    test, n = test_set(number=number)
    ind_positive = []
    ind_positive_uncertain = []
    ind_negative = []
    ind_hard = []
    ind_uncertain = []
    for i in tqdm(range(number)):
        f = test[i]
        X = IntervalVector(np.full((n[i],), -10), np.full((n[i],), 9))
        X_new, exist = algo(f, f.grad, X, **params)
        if (X_new != None) & exist:
            if X_new.w() < 1e-5:
                ind_positive += [i]
            else:
                ind_positive_uncertain += [i]
        if (X_new != None) & (not exist):
            if IntervalExtension(f)(X_new).M() < 1e-5:
                ind_uncertain += [i]
            else:
                ind_hard += [i]
        if (X_new == None):
            ind_negative += [i]
    return ind_positive, ind_positive_uncertain, ind_negative, ind_uncertain, ind_hard




def _krawczyk_operator_common_(f, f_ext, grad, slope_ext, dot_x, mid, X, slope=False):
    I = IntervalVector(np.identity(X.shape[0], dtype=object), np.identity(X.shape[0], dtype=object))
    
    A = grad(dot_x).astype(float)
    Y = np.linalg.pinv(A)

    if slope:
        S = slope_ext(X, mid)
    else:
        S = slope_ext(X)

    return (mid - Y @ f_ext(mid)) + ((I - Y @ S) @ (X - mid))

def krawczyk_method_common(f, grad, X_start, num_iters=10, tol=0, slope=None):    
    X = X_start
    i = 0

    f_ext = IntervalExtension(f)
    if slope != None:
        slope_ext = slope
        has_slope = True
    else:
        has_slope = False
        slope_ext = IntervalExtension(grad)
    
    X_old = X
    while (i < num_iters):
        mid = X.mid()
        dot_x = mid.a.astype(float)

        X_new = _krawczyk_operator_common_(f, f_ext, grad, slope_ext, dot_x, mid, X, has_slope)
        X = X & X_new

        if X == None:
            return None, False

        if (f_ext(X).M() <= tol):
            return X, True

        if X.is_eq(X_old):
            return X, False
        else:
            X_old = X
        i += 1
    return X, False

# def check_existence_krawczyk_common(f, grad, X_start, eps = dec.Decimal(1e-10), slope=None):    
#     X = X_start
#     X = IntervalVector(X.a - eps, X.b + eps)

#     f_ext = IntervalExtension(f)
#     if slope != None:
#         slope_ext = slope
#         has_slope = True
#     else:
#         has_slope = False
#         slope_ext = IntervalExtension(grad)
#     mid = X.mid()
#     dot_x = mid.a.astype(float)

#     X_new = _krawczyk_operator_common_(f, f_ext, grad, slope_ext, dot_x, mid, X, has_slope)
#     if (isin(X_new, X)):
#         return X, True
#     return X_start, False

def check_existence_krawczyk_common(f, grad, X_start, eps = dec.Decimal(1e-10), slope=None):    
    X = X_start
    X = IntervalVector(X.a - eps, X.b + eps)

    f_ext = IntervalExtension(f)
    if slope != None:
        slope_ext = slope
        has_slope = True
    else:
        has_slope = False
        slope_ext = IntervalExtension(grad)
    mid = X.mid()
    dot_x = mid.a.astype(float)

    
    A = grad(dot_x).astype(float)
    Y = np.linalg.pinv(A)
    Y, max_rank = is_max_rank_floor(Y)
    if max_rank:
        I = IntervalVector(np.identity(X.shape[0], dtype=object), np.identity(X.shape[0], dtype=object))
        if has_slope:
            S = slope_ext(X, mid)
        else:
            S = slope_ext(X)

        X_new = (mid - Y @ f_ext(mid)) + ((I - Y @ S) @ (X - mid))
        if (isin(X_new, X)):
            return X, True
    return X_start, False


def dissect_search_common(f, grad, X_start, steps=1000, algo=krawczyk_method_common, \
                          num_iters_algo=5, check_existence_algo=check_existence_krawczyk_common,\
                          tol=1e-4, slope=None, eps = dec.Decimal(1e-10)):
    l = [X_start]
    f_ext = IntervalExtension(f)
    to_dissect = False
    
    while (steps > 0) and (len(l) > 0):
        X = l.pop()

        if to_dissect:
            to_dissect = False
            
#             heur = f_ext(X)
#             heur = np.stack([np.abs(heur.a), np.abs(heur.b)])
#             h1 = heur.max(axis=0)
#             h2 = heur.min(axis=0)
#             ind = (h2 != 0)
#             heur = np.full(h1.shape, dec.Decimal(1))
#             heur[ind] = (h2[ind] / h1[ind])
#             ind = heur.argmax()   

            ind = X.coordinate_w().argmax()

#             l += X.dissect(ind, X.dot_mid()[ind])
#             X = l.pop()
#             if not contains_zero(f_ext(X)):
#                 X = l.pop()


            diss = X.dissect(ind, X.mid().a[ind])
            variant0 = f_ext(diss[0])
            variant1 = f_ext(diss[1])
            if not contains_zero(variant1):
                if contains_zero(variant0):
                    X = diss[0]
                else:
                    if len(l) > 0:
                        X = l.pop()
                    else:
                        return None, False
            else:
                if not contains_zero(variant0):
                    X = diss[1]
                else:
                    X = diss[1]
                    l += [diss[0]]

#                     variant0 = variant0[ind]
#                     variant1 = variant1[ind]
#                     heur = np.abs(np.array([variant0.a, variant0.b]))
#                     h0 = heur.min() / heur.max()
#                     heur = np.abs(np.array([variant1.a, variant1.b]))
#                     h1 = heur.min() / heur.max()
#                     if h0 > h1:
#                         X = diss[0]
#                         l += [diss[1]]
#                     else:
#                         X = diss[1]
#                         l += [diss[0]]

        tmp, is_tol = algo(f, grad, X, num_iters_algo, tol, slope)
        
        if is_tol:
            return check_existence_algo(f, grad, tmp, eps=eps, slope=slope)

        if tmp != None:
            l += [tmp]
            to_dissect = True
        steps -= 1
    if len(l) > 0:
        X = l.pop()
        return check_existence_algo(f, grad, X, eps=eps, slope=slope)
    else:
        return None, False

    
def dissect_algorithm_common(f, grad, X_start, steps=5000, num_iters=5, algo=krawczyk_method_common,\
                             check_existence_algo=check_existence_krawczyk_common,\
                             num_iters_algo=5, tol=1e-4, slope=None, eps=dec.Decimal(1e-10)):
    l = [X_start]
    certain = []
    uncertain = []
    start_steps = steps

    f_ext = IntervalExtension(f)
    to_dissect = False
    
    while (steps > 0) and (len(l) > 0):
        X = l.pop()
        if to_dissect:
            to_dissect = False

            ind = X.coordinate_w().argmax()


            diss = X.dissect(ind, X.mid().a[ind])
            variant0 = f_ext(diss[0])
            variant1 = f_ext(diss[1])
            if not contains_zero(variant1):
                if contains_zero(variant0):
                    X = diss[0]
                else:
                    if len(l) > 0:
                        X = l.pop()
                    else:
                        break
            else:
                if not contains_zero(variant0):
                    X = diss[1]
                else:
                    X = diss[1]
                    l += [diss[0]]

        tmp, is_tol = algo(f, grad, X, num_iters_algo, tol, slope)
        
        if is_tol:
            X_new, exist = check_existence_algo(f, grad, tmp, eps=eps, slope=slope)
            if exist:
                certain += [X_new]
            else:
                uncertain += [X_new]
        elif tmp != None:
            l += [tmp]
            to_dissect = True
        steps -= 1
    return certain, uncertain, l, start_steps - steps



def _gauss_seidel(X, dot_x, A, B, Y):
    n = len(X)
    diff_x = (X - dot_x).to_vec()

    A_cols = []
    exist = True
    for j in range(n):
        A_cols += [A[:, j]]

        
    X_new = []
    diff_x_new = []
    for i in range(n):
        Yi = IntervalVector.from_vector(Y[i, :])

        cum_sum = 0
        for j in range(i):
            cum_sum = cum_sum + (Yi @ A_cols[j]) * (diff_x_new[j])
        for j in range(i + 1, n):
            cum_sum = cum_sum + (Yi @ A_cols[j]) * (diff_x[j])

        xi = dot_x[i] - ((cum_sum - (Yi @ B)) / (Yi @ A_cols[i]))
        if not isin_1d(xi, X[i]):
            exist = False

        new = intersect_intervals(xi, X[i])

        if new == None:
            return None, False
        
        X_new += [new]
        diff_x_new += [new - dot_x[i]]
    return IntervalVector.from_vector(np.array(X_new)), exist

def _newton_step(f, f_ext, grad, grad_ext, X, dot_x, mid, slope=None):
    if slope is None:
        A = grad_ext(X)
    else:
        A = slope(X, mid)

    C = grad(dot_x).astype(float)
    Y = np.linalg.pinv(C)
    
    return _gauss_seidel(X, mid, A, -f_ext(mid), Y)

def newton_method(f, grad, X_start, num_iters=10, tol=1e-4, slope=None):
    X = X_start
    i = 0

    f_ext = IntervalExtension(f)
    grad_ext = IntervalExtension(grad)
    
    X_old = X
    while (i < num_iters):
        mid = X.mid()
        dot_x = mid.a

        X, exist = _newton_step(f, f_ext, grad, grad_ext, X, dot_x, mid, slope)

        if X == None:
            return None, False

        if (f_ext(X).M() <= tol):
            return X, True

        if X.is_eq(X_old):
            return X, False
        else:
            X_old = X
        i += 1
    return X, False

def check_existence_newton_common(f, grad, X_start, eps = dec.Decimal(1e-10), slope=None):    
    X = X_start
    X = IntervalVector(X.a - eps, X.b + eps)

    f_ext = IntervalExtension(f)
    if slope != None:
        slope_ext = slope
        has_slope = True
    else:
        has_slope = False
        slope_ext = IntervalExtension(grad)
    mid = X.mid()
    dot_x = mid.a.astype(float)
 
    if slope is None:
        A = grad_ext(X)
    else:
        A = slope(X, mid)

    C = grad(dot_x).astype(float)
    Y = np.linalg.pinv(C)
    Y, max_rank = is_max_rank_floor(Y)
    if max_rank:
        X_new, exist = _gauss_seidel(X, mid, A, -f_ext(mid), Y)
        if exist:
            return X_new, True
    return X_start, False

#     X_new, exist = _newton_step(f, f_ext, grad, slope_ext, X, dot_x, mid, slope=slope)
#     if exist:
#         return X_new, True
#     return X_start, False
