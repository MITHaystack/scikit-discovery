# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Soubhik Barari
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


"""Provides some example random walk functions in the expected format.

All random walk functions have the following signature:

(pos: tuple of input point, grid: bounds for walk, step_size: maximal step size) -> (output pos: tuple)

"""

import  random
from    numpy.random import multivariate_normal


def uniform_walk(pos, grid, step_size=None):
    '''
    A uniform random walk function
    
    @param pos: tuple of input point
    @param grid: bounds for walk
    @param step_size: maximal step size
    
    @return position tuple
    '''
    newpos = [float(random.randrange(grid[0]-grid[2]/2, grid[0]+grid[2]/2)), 
            float(random.randrange(grid[1]-grid[3]/2, grid[1]+grid[3]/2))]
            
    return newpos


def gaussian_walk(pos, grid, step_size=None):
    '''
    A gaussian random walk function
    
    @param pos: tuple of input point
    @param grid: bounds for walk
    @param step_size: maximal step size
    
    @return position tuple
    '''

    step = multivariate_normal(mean=[0,0], cov=[[grid[0],0],[0,grid[1]]], size=1)[0]
    newpos = [pos[0] + step[0], pos[1] + step[1]]       

    newpos = keep_in_bound(newpos, grid)
    return newpos


def keep_in_bound(pos, grid):
    '''
    Function for truncating and bounding the random walk to within the defined grid
    
    @param pos: tuple of the point to be checked
    @param grid: the bounds for limiting the walk
    
    @return position tuple after bounding the point
    '''
    pos[0] = max(pos[0], grid[0]-grid[2]/2)
    pos[0] = min(pos[0], grid[0]+grid[2]/2)
    pos[1] = max(pos[1], grid[1]-grid[3]/2)
    pos[1] = min(pos[1], grid[1]+grid[3]/2)
    return pos


