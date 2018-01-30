# The MIT License (MIT)
# Copyright (c) 2018 Massachusetts Institute of Technology
#
# Author: Guillaume Rongier
# We acknowledge support from NSF ACI1442997 (PI: V. Pankratius)
#                         and NASA AISTNNX15AG84G (PI: V. Pankratius)
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

#################################################################################
# Package imports
#################################################################################

import numpy as np
import math

from heapq import heappush, heappop

from numba import jit

#################################################################################
# Function definitions
#################################################################################

class PriorityQueue(object):

    def __init__(self, task_list, priority_list):
        count_list = list(range(len(task_list))) 
        # list of entries arranged in a heap
        self.pq_ = [list(i) for i in zip(priority_list, count_list, task_list)]
        # mapping of tasks to entries
        self.entry_finder_ = {task: entry for (task, entry) in zip(task_list, self.pq_)} 
        # unique sequence count
        self.counter_ = len(task_list) - 1
        
    def __repr__(self):
        return str(self.pq_)

    def add_task(self, task, priority = 0):
        'Add a new task or update the priority of an existing task'
        if task in self.entry_finder_:
            self.remove_task(task)
        self.counter_ += 1
        entry = [priority, self.counter_, task]
        self.entry_finder_[task] = entry
        heappush(self.pq_, entry)

    def remove_task(self, task):
        'Mark an existing task as REMOVED.  Raise KeyError if not found.'
        entry = self.entry_finder_.pop(task)
        entry[-1] = 'REMOVED'

    def pop_task(self):
        'Remove and return the lowest priority task. Raise KeyError if empty.'
        while self.pq_:
            priority, count, task = heappop(self.pq_)
            if task is not 'REMOVED':
                self.entry_finder_.pop(task, None)
                return task, priority
        raise KeyError('pop from an empty priority queue')
            
    def length(self):
        return len(self.entry_finder_)

    def is_empty(self):
        if len(self.entry_finder_) == 0:
            return True
        return False
    
    def empty(self):
        self.entry_finder_.clear()

@jit(nopython = True)
def get_four_neighborhood(j, 
                          i,
                          raster_height,
                          raster_width, 
                          gap = 1,
                          is_entire_planet_mapped = True):

    neighbors = [(-99999, -99999)]
    if j >= gap:
        neighbors.append((j - gap, i))
    elif is_entire_planet_mapped == True:
        neighbors.append((gap - j - 1, raster_width - 1 - i))
    if j < raster_height - gap:
        neighbors.append((j + gap, i))
    elif is_entire_planet_mapped == True:
        neighbors.append((2*raster_height - j - gap - 1, raster_width - 1 - i))
    if i >= gap:
        neighbors.append((j, i - gap))
    elif is_entire_planet_mapped == True:
        neighbors.append((j, raster_width - gap + i))
    if i < raster_width - gap:
        neighbors.append((j, i + gap))
    elif is_entire_planet_mapped == True:
        neighbors.append((j, gap - raster_width + i))
    
    return neighbors[1:]

@jit(nopython=True)
def haversine_distance_math(longitude_1, latitude_1, longitude_2, latitude_2, radius):
    rad_longitude_1 = longitude_1*math.pi/180.
    rad_latitude_1 = latitude_1*math.pi/180.
    rad_longitude_2 = longitude_2*math.pi/180.
    rad_latitude_2 = latitude_2*math.pi/180.
    
    diff_longitudes = rad_longitude_2 - rad_longitude_1
    diff_latitudes = rad_latitude_2 - rad_latitude_1
    
    root = math.sin(diff_latitudes/2.)**2 + math.cos(rad_latitude_1)*math.cos(rad_latitude_2)*math.sin(diff_longitudes/2.)**2
    return 2*radius*math.asin(math.sqrt(root))

@jit(nopython = True)
def get_quadratic_coefficients(current_cell_j, 
                               current_cell_i, 
                               time_array,
                               alive_cells,
                               velocity_array, 
                               longitude_array, 
                               latitude_array, 
                               planet_radius,
                               is_entire_planet_mapped = True):
    
    a = 0.
    b = 0.
    c = math.inf
    if velocity_array[current_cell_j, current_cell_i] != 0.:
        c = -1./(velocity_array[current_cell_j, current_cell_i]**2)
    
    raster_width = time_array.shape[2]
    raster_height = time_array.shape[1]
    neighboring_cells_1 = get_four_neighborhood(current_cell_j, current_cell_i, 
                                                raster_height, raster_width,
                                                is_entire_planet_mapped = is_entire_planet_mapped)
    neighboring_cells_2 = get_four_neighborhood(current_cell_j, current_cell_i, 
                                                raster_height, raster_width, 
                                                gap = 2, is_entire_planet_mapped = is_entire_planet_mapped)
    
    for neighboring_cell_1, neighboring_cell_2 in zip(neighboring_cells_1, neighboring_cells_2): 
    
        time_1 = time_array[0, neighboring_cell_1[0], neighboring_cell_1[1]]
        time_2 = time_array[0, neighboring_cell_2[0], neighboring_cell_2[1]]
        
        if (alive_cells[neighboring_cell_1] == 1 and
            alive_cells[neighboring_cell_2] == 1 and
            time_2 <= time_1):
            distance_1 = haversine_distance_math(longitude_array[neighboring_cell_1], 
                                                 latitude_array[neighboring_cell_1], 
                                                 longitude_array[current_cell_j, current_cell_i], 
                                                 latitude_array[current_cell_j, current_cell_i],  
                                                 planet_radius)
            distance_2 = haversine_distance_math(longitude_array[neighboring_cell_2], 
                                                 latitude_array[neighboring_cell_2], 
                                                 longitude_array[neighboring_cell_1], 
                                                 latitude_array[neighboring_cell_1],  
                                                 planet_radius)
            alpha = math.inf
            if distance_1 + distance_2 != 0.:
                alpha = 9./((distance_1 + distance_2)*(distance_1 + distance_2))
            d_time = (1./3.)*(4.*time_1 - time_2)
            a += alpha
            b -= 2.*alpha*d_time
            c += alpha*d_time*d_time
        elif alive_cells[neighboring_cell_1] == 1:
            distance = haversine_distance_math(longitude_array[neighboring_cell_1], 
                                               latitude_array[neighboring_cell_1], 
                                               longitude_array[current_cell_j, current_cell_i], 
                                               latitude_array[current_cell_j, current_cell_i],  
                                               planet_radius)
            alpha = math.inf
            if distance != 0.:
                alpha = 1/(distance*distance)
            a += alpha
            b -= 2.*alpha*time_1
            c += alpha*time_1*time_1

    return a, b, c

@jit(nopython = True)
def solve_quadratic_equation(a, b, c):
    delta = b*b - 4*a*c
    if delta >= 0. and a != 0.:
        # Fast marching calls for the largest possible solution to the quadratic equation
        return (-b + math.sqrt(delta))/(2*a)
    else:
        return math.inf #np.nan

@jit(nopython = True)
def compute_time(current_cell_j, 
                 current_cell_i, 
                 time_array,
                 alive_cells,
                 velocity_array, 
                 longitude_array, 
                 latitude_array, 
                 planet_radius,
                 is_entire_planet_mapped = True):
    a, b, c = get_quadratic_coefficients(current_cell_j, 
                                         current_cell_i, 
                                         time_array,
                                         alive_cells,
                                         velocity_array, 
                                         longitude_array, 
                                         latitude_array, 
                                         planet_radius,
                                         is_entire_planet_mapped)
    return solve_quadratic_equation(a, b, c)

def run_fast_marching(initiation_array,
                      velocity_array,
                      longitude_array,
                      latitude_array,
                      planet_radius,
                      stopping_time = None,
                      is_entire_planet_mapped = True,
                      turn_inf_to_nan = True):
    
    raster_width = initiation_array.shape[1]
    raster_height = initiation_array.shape[0]
    
    initial_cells = np.where((initiation_array == 0) &\
                             (~np.isnan(velocity_array)) &\
                             (velocity_array > 0.))
    
    time_array = np.full((3, raster_height, raster_width), np.nan)
    time_array[0,:,:] = np.inf
    time_array[0, (initiation_array == 0) &\
                  (~np.isnan(velocity_array)) &\
                  (velocity_array > 0.)] = 0.
    time_array[0, (initiation_array == 0) & (np.isnan(velocity_array))] = np.nan
    time_array[1:, (initiation_array == 0) &\
                   (~np.isnan(velocity_array)) &\
                   (velocity_array > 0.)] = np.array(initial_cells)

    if len(initial_cells[0]) > 0:
        
        initial_cells = list(zip(initial_cells[0], initial_cells[1]))
        initial_times = len(initial_cells)*[0]
        narrow_band = PriorityQueue(initial_cells, initial_times)

        alive_cells = np.zeros(initiation_array.shape, dtype = np.uint8)
        alive_cells[initiation_array == 0.] = 1

        while narrow_band.is_empty() == False:
            current_cell, current_time = narrow_band.pop_task()
            if stopping_time is not None and current_time > stopping_time:
                narrow_band.empty()
                time_array[0, alive_cells == 0] = np.inf
                time_array[1:, alive_cells == 0] = np.nan
            else:
                alive_cells[current_cell] = 1
                neighboring_cells = get_four_neighborhood(current_cell[0], current_cell[1], 
                                                          raster_height, raster_width,
                                                          is_entire_planet_mapped = is_entire_planet_mapped)
                for neighboring_cell in neighboring_cells:
                    if alive_cells[neighboring_cell] == 0:
                        neighbor_time = compute_time(neighboring_cell[0], 
                                                     neighboring_cell[1], 
                                                     time_array,
                                                     alive_cells,
                                                     velocity_array, 
                                                     longitude_array, 
                                                     latitude_array, 
                                                     planet_radius,
                                                     is_entire_planet_mapped)
                        if neighbor_time < time_array[0, neighboring_cell[0], neighboring_cell[1]]:
                            time_array[0, neighboring_cell[0], neighboring_cell[1]] = neighbor_time
                            time_array[1:,
                                       neighboring_cell[0],
                                       neighboring_cell[1]] = time_array[1:,
                                                                         current_cell[0],
                                                                         current_cell[1]]
                            narrow_band.add_task(neighboring_cell, neighbor_time)
    if turn_inf_to_nan == True:
        time_array[np.isinf(time_array)] = np.nan

    return time_array