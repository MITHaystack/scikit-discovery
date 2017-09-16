# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Cody Rude
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

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
from scipy.optimize import root
from skdaccess.framework.data_class import DataFetcherBase, TableWrapper

from skdiscovery.utilities.patterns.astro_tools import nfw, move_point


class CatalogGenerator(DataFetcherBase):
    '''
    *In Development* Generates galaxy catalogs for use in DiscoveryPipeline
    '''

    def __init__(self, ap_paramList, ra1,dec1, ra2, dec2, background_density, z):
        '''
        @param ap_paramList[seed]: Seed for random number generator
        @param ra1:  Left right ascension
        @param dec1: Bottom declination
        @param ra2:  Right right ascension
        @param dec2: Top declination
        @param background_density: galaxy background density in galaxies/square degree
        @param z: Redshift of galaxy cluster
        '''

        self.ra1 = ra1
        self.dec1 = dec1
        self.ra2 = ra2
        self.dec2 = dec2
        self.background_density = background_density
        self.z = z

        self.__H0 = 70
        self.__Omega_m = 0.3

        self.__R0 = 2
        self.__Rs = 0.15 / 0.7
        self.__Rcore = 0.1 / 0.7
        self.__norm = 1/quad(nfw,0,self.__R0,args=(1,self.__Rs,self.__Rcore))[0]

        super(CatalogGenerator, self).__init__(ap_paramList)
    

    def output(self):
        ''' 
        Generates galaxy catalog

        @return DataWrapper: Table data wrapper of galaxy catalog
        '''

        if len(self.ap_paramList) > 0:
            seed = self.ap_paramList[0]()
            np.random.seed(seed)

        else:
            seed = None

        # Generating a background for the whole sky, and then selecting galaxies
        # from that background
        num_background_galaxies = round(self.background_density * (4 * np.pi) * (180./np.pi)**2)

        # Generate background galaxies over an entire sphere
        full_ra = np.random.rand(num_background_galaxies) * 360.
        full_dec = np.arccos(2*np.random.rand(num_background_galaxies)-1) * 180. / np.pi - 90.

        # Make data frame and select galaxies in our box
        galaxies = pd.DataFrame.from_dict({'RA' : full_ra, 'Dec' : full_dec})

        galaxies = galaxies[ np.logical_and.reduce( (galaxies['RA'] > self.ra1, galaxies['RA'] < self.ra2,
                                                     galaxies['Dec'] > self.dec1, galaxies['Dec'] < self.dec2) ) ]

        galaxies['Cluster_ID'] = -1        


        # Now generating a cluster
        num_galaxies = 100

        cluster_ra = np.mean([self.ra1,self.ra2])
        cluster_dec = np.mean([self.dec1, self.dec2])

        cosmo = FlatLambdaCDM(self.__H0, self.__Omega_m)

        distance = cosmo.comoving_distance(self.z).to('Mpc').value

        radial_positions = self.inverse_nfw_cumulative(np.random.rand(num_galaxies))
        bearings = np.random.rand(num_galaxies) * 360.

        # Converting radial comoving distances to angular seperations
        angular_distances = radial_positions / distance * 180./np.pi

        # Generate the RA and Dec for each galaxy by moving it in random direction
        # the radial distance away from the center
        
        cluster_ra, cluster_dec = move_point(cluster_ra, cluster_dec, angular_distances, bearings)
        cluster_galaxies = pd.DataFrame.from_dict( {'RA' : cluster_ra, 'Dec' : cluster_dec} )

        cluster_galaxies['Cluster_ID'] = int(0)

        galaxies = pd.concat([galaxies, cluster_galaxies], axis=0).set_index(np.arange(len(galaxies) + len(cluster_galaxies)))

        data_wrapper = TableWrapper({'Cluster_Catalog_01' : galaxies}, default_columns=['RA','Dec'])
        
        return data_wrapper




    def nfw_cumulative(self,R):
        ''' 
        Cumulative radial NFW distribution

        @param R: Radius
        @return float: Probability of being within R
        '''
        
        R0 = self.__R0
        norm = self.__norm
        Rs = self.__Rs
        Rcore = self.__Rcore

        def integrate(z):
            return quad(nfw,0,z,args=(norm,Rs,Rcore))[0]

        if np.isscalar(R):
            R = np.array([float(R)])
        else:
            R = R.astype(np.float)

        res = np.piecewise(R, [R < 0.0, np.logical_and(R >= 0.0, R < R0), R >= R0], 
                           [lambda z: z, 
                            lambda z: np.fromiter(map(integrate,z),np.float),
                            lambda z: z/R0])
        if len(res) == 1:
            return res[0]
        else:
            return res

    def inverse_nfw_cumulative(self, p):
        '''
        inverse of radial nfw cumulative distribution

        @param p: Probability
        @return float: Radius corresponding to probability p 
        '''
        R0 = self.__R0

        def get_root(p):
            return root(lambda z: self.nfw_cumulative(z) - p,R0/2).x[0]
        if np.isscalar(p):
            get_root(p)
        else:
            return np.fromiter(map(get_root,p),np.float, count=len(p))
