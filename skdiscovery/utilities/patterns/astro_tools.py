# The following license applies to angular_seperation, abs_mag,
# app_mag, and move_point

# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Victor Pankratius, Justin Li, Cody Rude
# This software is part of the NSF DIBBS Project "An Infrastructure for
# Computer Aided Discovery in Geoscience" (PI: V. Pankratius) and 
# NASA AIST Project "Computer-Aided Discovery of Earth Surface 
# Deformation Phenomena" (PI: V. Pankratius)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# The functions z_to_v, v_to_z, lf, dlf, cdf_dlf, and inv_cdf_dlf are
# adapted from R code written by Cody Rude. The function NFW was
# adapted from python code written by Cody Rude. The original versions
# of these are licensed under the MIT license below. Any changes made
# to the original version are copyright MIT and follow the license and
# copyright listed above.


# The MIT License (MIT)
# Copyright (c) 2016 Cody Rude
#
# Author: Cody Rude
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import numpy as np
import math
from astropy.constants import c
from astropy import units as u
from astropy.coordinates import SkyCoord as SC
from astropy.cosmology import FlatLambdaCDM

from scipy.optimize import root
from scipy.integrate import quad


def z_to_v(z):
    '''
    Convert redshift to km/s assuming shift is due to velocity using special relativity

    @param z: Redshift

    @return speed in km/s assuming shift is due to motion using special relativity
    '''
    ckms = c.to('km/s').value
    return (ckms * ( (z+1)**2 - 1) / ( (z+1)**2 + 1 ))


def v_to_z(v):
    '''
    Convert km/s to redshift assuming all are using special relativity

    @param v: velocity in km/s

    @return Redshift of object with speed in km/s
    '''
    ckms = c.to('km/s').value    
    return np.sqrt((1 + v/ckms) / (1 - v/ckms)) - 1


def angular_separation(ra1, dec1, ra2, dec2):
    '''
    Angular seperation between two objects via the haversine formula. 
    
    All inputs are in degrees.

    Formula obtained from http://www.movable-type.co.uk/scripts/gis-faq-5.1.html

    Formula originally presented in
    R.W. Sinnott, "Virtues of the Haversine", Sky and Telescope, vol. 68, no. 2, 1984, p. 159

    @param ra1: Right Ascention of first object (degrees)
    @param dec1: Declination of first object (degrees)
    @param ra2: Right Ascention of second object (degrees)
    @param dec2: Declination of second object (degrees)

    @return angular seperation between two objects
    '''

    ra1 = ra1 * np.pi / 180.
    dec1 = dec1 * np.pi / 180.
    ra2 = ra2 * np.pi / 180.
    dec2 = dec2 * np.pi / 180.

    dra  = ra2 - ra1
    ddec = dec2 - dec1

    a = np.sin(ddec/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dra/2)**2
    return 2*np.arcsin(np.minimum(1,np.sqrt(a))) * 180. / np.pi
    


def move_point(ra,dec,ang_dist,bearing):
    '''
    Move a point along a great circle at a particular bearing

    All inputs are in degrees
    The formula was obtained from http://www.movable-type.co.uk/scripts/latlong.html

    @param ra: Starting right ascension
    @param dec: Starting declination
    @param ang_dist: Angular distance to travel
    @param bearing: Direction to travel (0 is north, 90 is positive RA)

    @return tuple containing updated ra and dec
    '''
    
    ra = ra * np.pi / 180.
    dec = dec * np.pi / 180.
    ang_dist = ang_dist * np.pi / 180.
    bearing = bearing * np.pi / 180
    
    new_dec = np.arcsin( np.sin(dec) * np.cos(ang_dist) + np.cos(dec)*np.sin(ang_dist)*np.cos(bearing))
    new_ra = ra + np.arctan2(np.sin(bearing) * np.sin(ang_dist) * np.cos(dec), 
                             np.cos(ang_dist) - np.sin(dec)*np.sin(new_dec))
    
    new_dec *= 180. / np.pi
    new_ra *= 180. / np.pi
    
    new_ra %= 360.
    
    return new_ra, new_dec
    

def abs_mag(app_mag, z):
    ''' 
    Get the absolute magnitude from apparent magnitude
    
    Assumes concordance cosmology. No kcorrection is applied.

    @param app_mag: Apparent magnitude
    @param z: Redshift

    @return absolute magnitude of object at z
    '''
    
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    distmod = 5*np.log10(cosmo.luminosity_distance(z).to('pc').value)-5
    return app_mag - distmod


def app_mag(abs_mag, z):
    ''' 
    Get the apparent magnitude from absolute magnitude.

    Assumes concordance cosmology. No kcorrection is assumed.

    @param abs_mag: Absolute magnitude
    @param z: Redshift

    @return apparent magnitude of object at z
    '''

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    distmod = 5*np.log10(cosmo.luminosity_distance(z).to('pc').value)-5
    return abs_mag + distmod


def nfw(R, norm_constant, Rs, Rcore):
    '''
    2D Navarro-Frenk-White surface radial profile probability density

    See:
    Navarro, J. F., Frenk, C. S., & White, S. D. M. 1996, ApJ, 462, 563
    Bartelmann, M., A&A, 1996, 313, 697
    Rykoff, E.S., et al., ApJ, 746, 178

    @param R: Radius
    @param norm_constant: Normalization constant
    @param Rs: Scale radius
    @param Rcore: Since NFW profile diverges at R=0, the value at the center is held fixed starting at Rcore

    @return probability density of profile at R
    '''

    def compute_nfw(R):
        if R < Rcore:
            R2 = Rcore
        else:
            R2 = R

        if (R2/Rs) < 0.999:
            func = (1 - 2 * math.atanh(math.sqrt((1 - R2/Rs) / (R2/Rs + 1))) / math.sqrt(1 - (R2/Rs)**2)) / ((R2/Rs)**2 - 1)

        # There are some computational issues as R2/Rs -> 1, using taylor expansion of the function
        # around this point
        elif (R2/Rs) < 1.001:
            func = -(20/63)*((R2/Rs)**3 - 1) + (13/35)*((R2/Rs)**2 - 1) - (2/5)*(R2/Rs) + (11/15)

        else:
            func = (1 - 2 * math.atan(math.sqrt((R2/Rs - 1) / (R2/Rs + 1))) / math.sqrt((R2/Rs)**2 - 1)) / ((R2/Rs)**2 - 1)


        return norm_constant * 2 * math.pi * R * func

    if np.isscalar(R):
        return compute_nfw(R)

    else:
        return np.fromiter(map(compute_nfw, R), np.float, count = len(R))


def lf(x, A, mstar, alpha):
    ''' 
    Schechter function

    @param x: magnitude
    @param A: Scale factor
    @param mstar: Knee of distribution
    @param alpha: Faint-end turnover
    
    @return float: Schecter function at magnitude x
    '''
    return(A * 10**(-0.4 * (x - mstar) * (alpha + 1)) * np.exp(-10**(-0.4 * (x - mstar)) ))


def dlf(x, A, m1, a1, m2, a2):
    ''' 
    double Schechter function. Second LF is set to be 2*A of first LF.

    @param x: magnitude
    @param A: Scale factor
    @param m1: Knee of distribution 1
    @param a1: Faint-end turnover of first lf
    @param m2: Knee of distribution 2
    @param a2: Faint-end turnover of second lf
    
    @return float: Double Schecter function at magnitude x
    '''
    
    return(lf(x, A, m1, a1) + lf(x, 2*A, m2, a2))


def cdf_dlf(x, A, m1, a1, m2, a2, start=-26):
    ''' 
    Cumulative  Schechter function. Second LF is set to be 2*A of first LF.

    @param x: magnitude
    @param A: Scale factor
    @param m1: Knee of distribution 1
    @param a1: Faint-end turnover of first lf
    @param m2: Knee of distribution 2
    @param a2: Faint-end turnover of second lf
    @param start: Brightest magnitude
    
    @return Probability that galaxy has a magnitude greater than x
    '''
    def integrate(in_x):
        return quad(dlf, start,in_x,args=(A,m1,a1,m2,a2))[0]

    if np.isscalar(x):
        x = np.array([x])

    return np.fromiter(map(integrate,x),np.float,count=len(x))


def inv_cdf_dlf(p, A, m1, a1, m2, a2, start=-26, end=-15):

    ''' 
    Inverse Cumulative Schechter function. Second LF is set to be 2*A of first LF.

    @param p: probability
    @param A: Scale factor
    @param m1: Knee of distribution 1
    @param a1: Faint-end turnover of first lf
    @param m2: Knee of distribution 2
    @param a2: Faint-end turnover of second lf
    @param start: Brightest magnitude
    @param end: Faintest possible magnitude
    
    @return Magnitude associated with cdf probability p
    '''
    def get_root(p):
        return root(lambda x: cdf_dlf(x,A,m1,a1,m2,a2,start)-p, (start + end)/2).x[0]

    if np.isscalar(p):
        return get_root(p)
    else:
        return np.fromiter(map(get_root,p),np.float,count=len(p))
