#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# This file is part of the PyNCulture project, which aims at providing tools to
# easily generate complex neuronal cultures.
# Copyright (C) 2017 SENeC Initiative
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# The tools present in this file are modified from the homcoord.py file by
# John W. Shipman (see __credits__)
#
# Contributor: S. Bottani

"""
2D Homogeneous coordinates with transformations

With addition of the SEGMENT class
S. Bottani
"""

__credits__ = ['http://www.nmt.edu/tcc/help/lang/python/examples/homcoord/']


import sys
from math import *

import numpy as np


# --------- #
# Constants #
# --------- #

RAD_45  = pi / 4.  # 45 degrees in radians
RAD_90  = pi / 2.  # 90 degrees in radians
RAD_180 = pi       # 180 degrees in radians
TWO_PI  = 2. * pi  # 360 degrees in radians


# ----------- #
# Point class #
# ----------- #

class Pt(object):

    '''
    Represents a homogeneous coordinate in 2-space.

    Exports:
    Pt(*coords):
      [ coords is a 2-sequence or a 1-sequence containing a
        2-sequence ->
          return a new Pt instance representing those two
          values as x and y, respectively
        coords is a 1-sequence containing a 3-sequence ->
          return a new Pt instance representing those values
          as x, y, and w, respectively ]
    .xy:
      [ return a 2-tuple with the homogenized x and y values ]
    .x:    [ return the homogenized x coordinate ]
    .y:    [ return the homogenized y coordinate ]
    .dist(other):
      [ other is a Pt instance ->
          return the distance between self and other ]
    .bearing(p):
      [ p is a Pt instance ->
          return the Cartesian angle in radians from self to p ]
    .radial(d, bearing):
      [ (d is a distance) and (bearing is an angle in radians) ->
          return the location at that distance and bearing as
          a Pt instance ]
    .toPolar():
      [ return self in polar coordinates as a Polar instance ]
    .__str__():  [ return self as a string ]
    .__add__(self, other):
      [ other is a Pt instance ->
          return a new Pt instance whose coordinates are the
          sum of self's and other's ]
    .__sub__(self, other):
      [ other is a Pt instance ->
          return a new Pt instance whose coordinates are the
          self's minus other's ]
    .__cmp__(self, other):
      [ if self and other are the same point ->
          return 0
        else -> return a nonzero value ]
    State/Invariants:
    .v     [ a numpy 3-element vector [x, y, W] ]
    '''

    def __init__ ( self, *args ):
        '''Constructor.
        '''
        #-- 1 --
        # [ if args has one element containing exactly two values ->
        #     x  :=  args[0][0]
        #     y  :=  args[0][1]
        #     w  :=  1.0
        #   else if args has one element containing three values ->
        #     x  :=  args[0][0]
        #     y  :=  args[0][1]
        #     w  :=  args[0][2]
        #   else if args has exactly two values ->
        #     x  :=  args[0]
        #     y  :=  args[1]
        #     w  :=  1.0
        #   else -> raise ValueError ]
        if len(args) == 1:
            value = args[0]
            if len(value) == 2:
                x, y = value
                w = 1.0
            else:
                x, y, w = value
        else:
            x, y = args
            w = 1.0
        #-- 2 --
        # [ self.v  :=  a 3-element numpy vector (x, y, w) as
        #                 type float
        #   self.__inverse  :=  None ]
        w=float(w) #to force floating division in x,y properties
        self.v = (x, y, w)

    @property
    def xy(self):
        '''
        Return (x,y)
        '''
        w = self.v[2]
        return (self.v[0]/w, self.v[1]/w)

    @property
    def x(self):
        '''
        Return the abscissa.
        '''
        return self.v[0]/self.v[2]

    @property
    def y(self):
        '''Return the ordinate.
        '''
        return self.v[1]/self.v[2]

    def apply(self,f):
        """:return: Pt obtained by appying function f to x and y"""
        return Pt(f(self.x),f(self.y))

    def dist(self, other):
        '''Return the distance between self and other.
        '''
        (dx,dy) = (self-other).xy
        return sqrt ( dx*dx + dy*dy )

    def bearing(self, p):
        '''
        What is the bearing angle from self to p?
        '''
        return atan2(p.y-self.y, p.x-self.x)

    def radial(self, d, bearing):
        '''
        Return the point at a given distance and bearing.
        '''
        return Pt(self.x + d*cos(bearing),
                  self.y + d*sin(bearing) )

    def toPolar(self):
        '''Convert to polar coordinates.
        '''
        x, y = self.xy
        return Polar(sqrt(x*x + y*y), atan2(y, x))

    def __str__(self):
        '''Return a string representation of self.
        '''
        return "(%.4g, %.4g)" % self.xy

    def __repr__(self):
        return str(self)

    def __add__(self, other):
        '''Add two points.
        '''
        return Pt(self.x+other.x, self.y+other.y)

    def __sub__(self, other):
        '''Subtract two points.
        '''
        return Pt(self.x-other.x, self.y-other.y)

    def __mul__(self, scale):
        '''Multiply by scalar.
        '''
        return Pt((self.x, self.y, 1./scale))

    def __div__(self, scale):
        '''Multiply by scalar.
        '''
        return Pt((self.x, self.y, scale))

    def __cmp__(self, other):
        '''Compare two points.
        '''
        return cmp(self.xy, other.xy)


# -------------------- #
# Transformation class #
# -------------------- #

class Xform(object):

    '''
    Represents an arbitrary homogeneous coordinate transform.

    Exports:
    Xform(m):
      [ m is a 3x3 transform matrix as a array, or
        a sequence that array() will accept as a 3x3
        array ->
          return a new Xform instance representing that
          transform ]
    .apply(p):
      [ p is a Pt instance ->
          return a new Pt instance representing p transformed
          by self ]
    .invert(p):
      [ p is a Pt instance ->
          return a new Pt instance pp such that
          self.apply(pp) == p ]
    .inverse():
      [ return the inverse of self as an Xform instance ]
    .compose(t):
      [ t is an Xform instance ->
          return a new Xform representing the composition of
          self followed by t ]
    .offset():
      [ return the net offset that self will shift the origin,
        as a Pt instance ]
    .angle():
      [ return the net angle that self will rotate the unit
        vector from (0,0) to (1,1) ]
    .mag():
      [ return the net magnification that self will apply to the
        unit vector ]
    .__str__(self):
      [ return a string representation of self ]
    State/Invariants:
    self._m:
      [ a 3x3 array representing the argument passed
        to the constructor ]
    self._mInverse:
      [ the inverse of self._m or None ]
    self.__offset:
      [ the net translation of self or None ]
    self.__angle:
      [ the net rotation of self or None ]
    self._mag:
      [ the net uniform scaling of self or None ]
    ORIGIN:      [ the origin as a Pt instance ]
    UNIT:        [ a point 1.0 along the line x=y ]
    '''

    ORIGIN = Pt(0,0)
    UNIT = ORIGIN.radial(1.0, RAD_45)

    def __init__ ( self, m ):
        '''Constructor.
        '''
        #-- 1 --
        # [ if the type of m is ndarray ->
        #     self._m  :=  m
        #   else if m is acceptable as an argument to
        #   array() ->
        #     self._m  :=  array(m)
        #   else -> raise Exception ]
        self._m = m
        #-- 2 --
        self._mInverse = None

    def apply ( self, p ):
        '''Transform a point.
        '''
        #-- 1 --
        # [ pp  :=  a array representing the dot product of
        #           self._m and p.v ]
        pp = dot(self._m, p.v)

        #-- 2 --
        # [ return a Pt instance representing pp.v ]
        return Pt(pp)

    def __call__(self, p):
        return self.apply(p)

    def invert ( self, p ):
        '''Return p transformed by the inverse of self, as a Pt.
        '''
        return self.inverse().apply ( p )

    def inverse ( self ):
        '''Return the inverse transform as an Xform.
        '''
        #-- 1 --
        # [ if self._mInverse is None ->
        #     self._mInverse  :=  matrix inverse of self._m
        #   else -> I ]
        if self._mInverse is None:
            self._mInverse = linalg.inv ( self._m )

        #-- 2 --
        return Xform(self._mInverse)

    def __str__(self):
        '''Display self as a string
        '''
        #-- 1 --
        return ( "<Xform(xlate(%s), rotate(%.1fdeg), "
                 "mag(%.1f)>" %
                 (self.offset(), degrees(self.angle()), self.mag()) )

    def compose ( self, t2 ):
        '''Return the composition of two transforms.
        '''
        return Xform ( dot ( t2._m, self._m ) )

    def __mul__ ( self, other ):
        '''Implement '*'
        '''
        return self.compose(other)

    def offset(self):
        return self(self.ORIGIN )

    def angle(self,angle=RAD_45):
        """
        :param angle: angle in radians of a unit vector starting at origin
        :return: float bearing in radians of the transformed vector
        """
        pt=Polar(1.0,angle).toCartesian()
        pt=self(pt)-self.offset()
        return atan2(pt.y,pt.x)


    def mag(self):
        '''Return the net (uniform) scaling of this transform.
        '''
        return self(self.ORIGIN ).dist(self(self.UNIT))


# ----------- #
# Polar class #
# ----------- #

class Polar(object):
    '''
    Represents a point in polar coordinates.

    Exports:
    Polar(r, theta):
      [ r and theta are numbers ->
          return a new Polar instance representing radius r
          and angle theta ]
    .r, .theta:  [ as passed to constructor ]
    .toCartesian():
      [ return self in Cartesian coordinates as a Pt instance ]
    .__str__():
      [ return self as a string "(r, theta)" ]
    '''

    def __init__(self, *p):
        '''Constructor
        '''
        self.r, self.theta = argPair(*p)

    def toCartesian(self):
        '''Return self in rectangular coordinates as a Pt instance.
        '''
        return Pt(self.r * cos(self.theta),
                  self.r * sin(self.theta))

    def __str__(self):
        '''Return self as a string.
        '''
        return ( "(%.4g, %.4gd)" %
                 (self.r, degrees(self.theta)) )


# ---------- #
# Line class #
# ---------- #

class Line(object):

    '''
    Represents a geometric line.

    Exports:
    Line(a, b, c):
      [ a, b, and c are floats ->
          return a Line instance representing ax+by+c=0 ]
    .a, .b, .c:  [ as passed to constructor, read-only ]
    .__str__(self):   [ return self as a string ]
    .intersect(other):
      [ other is a Line instance ->
          if self and other intersect ->
            return the intersection as a Pt
          else -> raise ValueError ]
    Line.twoPoint(p1, p2):       # Static method
      [ p1 and p2 are Pt instances ->
          if p1 and p2 are distinct ->
            return a Line instance representing the line that
            intersects p1 and p2
          else -> raise ValueError ]
    Line.pointBearing(p, bears):   # Static method
      [ (p is a Pt instance) and
        (bears is a Cartesian bearing in radians) ->
          return the line through p at bearing (bears) ]
    '''

    def __init__(self, a, b, c):
        '''Constructor.
        '''
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)

    def __str__(self):
        '''Return a string representing self.
        '''
        return "%.4gx + %.4gy + %.4g = 0" % (self.a, self.b, self.c)

    def intersect(self, other):
        '''Where do lines self and other intersect?
        '''
        #-- 1 --
        # [ if self and other have the same slope ->
        #     raise ValueError
        #   else -> I ]
        if self.a * other.b == other.a * self.b:
            raise ValueError("Lines have the same slope.")
        #-- 2 --
        # [ x, y  :=  solution to the simultaneous linear equations
        #       (self.a * x + self.b * y = -self.c) and
        #       (other.a * x + other.b * y = -other.c) ]
        a = array ( ( (self.a, self.b), (other.a, other.b) ) )
        b = array ( (-self.c, -other.c) )
        x, y = linalg.solve(a,b)

        #-- 3 --
        return Pt(x, y)

    @staticmethod
    def twoPoint(p1, p2):
        '''Find the equation of a line between two points.
        '''
        #-- 1 --
        # [ if p1 and p2 coincide ->
        #     raise ValueError
        #   else ->
        #     x1  :=  abscissa of p1
        #     y1  :=  ordinate of p1
        #     x2  :=  abscissa of p2
        #     y2  :=  ordinate of p2 ]
        if p1 == p2:
            raise ValueError("Points are not distinct.")
        else:
            x1, y1 = p1.xy
            x2, y2 = p2.xy
        #-- 2 --
        # [ if x1 == x2 ->
        #     return a vertical line through x1
        #   else ->
        #     m  :=  (y2-y1)/(x2-x1) ]
        if x1 == x2:
            return Line(1.0, 0.0, -x1)
        else:
            m = (y2-y1)/(x2-x1)
        #-- 4 --
        # [ return a new Line instance having a=(-m), b=1, and
        #   c=(m*x1-y1) ]
        return Line(-m, 1.0, (m*x1-y1))

    @staticmethod
    def pointBearing(p, bears):
        '''Line through p at angle (bears).
        '''
        #-- 1 --
        # [ angle  :=  angle normalized to [0,180) degrees
        #   px  :=  abscissa of p
        #   py  :=  ordinate of p ]
        angle = bears % RAD_180
        px, py = p.xy

        #-- 2 --
        # [ if angle == RAD_90 ->
        #     return a Line with a=1.0, b=0.0, and c=-p.x
        #   else ->
        #     m  :=  tan(angle) ]
        if angle == RAD_90:
            return Line(1.0, 0.0, -px)
        else:
            m = tan(angle)
        #-- 3 --
        # [ return a Line with a=m, b=-1.0, and c=(-m*px + py) ]
        return Line(m, -1.0, py - m*px)

def argPair(*p):
    '''Process a pair of values passed in various ways.

      [ if len(p) is 2 ->
            return (p[0], p[1])
        else if p is a single non-iterable ->
          return (p[0], p[0])
        else if p is an iterable with two values ->
            return (p[0][0], p[0][1])
        else if p is an iterable with one value ->
            return (p[0][0], p[0][0])
    '''
    #-- 1 --
    if len(p) == 2:
        return (p[0], p[1])

    #-- 2 --
    it = p[0]
    if not hasattr(it, "__iter__"):
        return(it, it)

    #-- 3 --
    # [ p is an iterable ->
    #     values  :=  all values from p[0] ]
    values = [ x
               for x in p[0] ]

    #-- 4 --
    if len(values) == 1:
        return (values[0], values[0])
    return (values[0], values[1])


# ------------- #
# Segment class #
# ------------- #

class Segment(object):
    '''Represents a geometric Segment.

      Exports:
        Segment(a, b):
          [ a, b, are Points ->
              returns a Segment instance representing the segment between points a and b
              .length is the segment length
              .dir the unit vector in the direction a to b]

        .__str__(self):   [ return self as a string ]
        .intersect(other):
          [ other is a Segment instance ->
              if self and other intersect ->
                return the intersection as a Pt
              else -> raise ValueError ]

    '''

    def __init__(self, a, b):
        '''Constructor.
        '''
        self.a = a # a, and b Pt
        self.b = b
        self.length=np.sqrt((b.x-a.x)**2+(b.y-a.y)**2)
        self.dir=np.array((b.x-a.x),(b.y-a.y))/self.length

    def __str__(self):
        '''Return a string representing self.
        '''
        return "(%.4gx , %.4gy )" % (self.a, self.b)

    def intersect(self, other):
        ''' If segments self and other intersect return the intersection point as
        a Pt instance
        '''
        #-- 1 --
        # [ if self and other have the same slope ->
        #     raise ValueError
        #   else -> I ]
        if self.dir == other.dir:
            return
        #    raise ValueError("Lines have the same slope.")
        #-- 2 --
        # [ x, y  :=  solution to the simultaneous linear equations
        #       (self.a * x + self.b * y = -self.c) and
        #       (other.a * x + other.b * y = -other.c) ]


        line1=Line.twoPoints(self.a,self.b)
        line2=Line.twoPoints(other.a,other.b)

        intersection=line1.intersect(line2) # returns 0 it the 2 lines do not intersect

        # Check if intersection is in the segments

        self_in_segment = is_in_segment(intersection, self)
        other_in_segment = is_in_segment(intersection,other)
        if self_in_segment and other_in_segment:
            return intersection
        return False


# ----- #
# Tools #
# ----- #

def dot(a,b):
    '''
    Generalized dot product.
    '''
    try: #vector*vector
        return sum(map( operator.mul, a, b))
    except:
        pass
    try: #matrix*vector
        return [dot(line,b) for line in a]
    except:
        pass
    #matrix*matrix
    res=[dot(a,col) for col in zip(*b)]
    return map(list, zip(*res))


def normAngle(theta):
    '''
    Normalize an angle in radians to [0, 2*pi).
    '''
    return theta % TWO_PI


def is_in_segment(pt,seg):
    '''
    Returns True if a point is in the segment.
    '''

    scaling_x=np.abs((pt.x-seg.b.x)/(pt.x-seg.a.x))
    scaling_y=np.abs((pt.y-seg.b.y)/(pt.y-seg.a.y))

    if scaling_x != scaling_y :
        return False

    if scaling_x >1. :
        return False

    return True


# TRANSFORM FUNCTIONS

def Xlate(*p):
    '''
    Create a translation transform.
    '''
    #-- 1 --
    # [ dx  :=  first value from p
    #   dy  :=  second value from p ]
    dx, dy = argPair ( *p )

    #-- 2 --
    return Xform ( [ (1, 0, dx),
                     (0, 1, dy),
                     (0, 0, 1)  ] )


def Xscale(*p):
    '''
    Create a scaling transform.
    '''
    #-- 1 --
    # [ if p is a single value or single-valued iterable ->
    #     sx  :=  that value
    #     sy  :=  that value
    #   else ->
    #     sx  :=  the first value from p
    #     sy  :=  the second value from p ]
    sx, sy = argPair ( *p )

    #-- 2 --
    # [ return an Xform for scaling x by sx and scaling y by sy ]
    return Xform ( [ (sx, 0,  0),
                     (0,  sy, 0),
                     (0,  0,  1) ] )


def Xrotate(theta):
    '''Create a rotation transform.
    '''
    #-- 1 --
    sint = sin(theta)
    cost = cos(theta)

    #-- 2 --
    return Xform ( [ (cost, -sint, 0),
                     (sint, cost,  0),
                     (0,    0,     1) ] )


def Xrotaround ( p, theta ):
    '''Rotation of theta radians around point p.
    '''
    #-- 1 --
    # [ t1  :=  an Xform that translates point p to the origin
    #   r  :=  an Xform that rotates theta radians around the origin
    #   t2  :=  an Xform that translates the origin to point p ]
    t1 = Xlate ( [ -v
                   for v in p.xy ] )
    r = Xrotate ( theta )
    t2 = Xlate ( p.xy )

    #-- 2 --
    # [ return an Xform instance representing t1, then r, then t2 ]
    return t1.compose(r).compose(t2)
