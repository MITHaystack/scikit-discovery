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

################################################################################
# Package imports
################################################################################

import csv

import math
import numpy as np

from collections import defaultdict
from collections import deque
from ast import literal_eval as make_tuple

from skdiscovery.utilities.planetary.geographical_computation import haversine_distance_math

################################################################################
# Function definitions
################################################################################

# Computation of the neighbors and threshold targets

def get_target_types_at_cells(target_arrays):
    target_types_at_cells = defaultdict(set)

    for index, array in enumerate(target_arrays):
        cells_to_visit = list(zip(*np.where(array[0,:,:] == 0)))
        for cell in cells_to_visit:
            target_types_at_cells[cell].add(index)
            
    return target_types_at_cells

def identify_neighbors(cells, target_types_at_cells, target_arrays, time_limit = math.inf):
    neighbors = defaultdict(list)

    for target in cells:
        target_neighbors = defaultdict(list)
        for k, k_array in enumerate(target_arrays):
            if (k not in target_types_at_cells[target]
                and math.isnan(k_array[0, target[0], target[1]]) == False):
                neighbor = (int(k_array[1, target[0], target[1]]),
                            int(k_array[2, target[0], target[1]]))
                target_neighbors[neighbor].append(k_array[0, target[0], target[1]])

                for l in target_types_at_cells[target]:
                    if (target_arrays[l][1, neighbor[0], neighbor[1]],
                        target_arrays[l][2, neighbor[0], neighbor[1]]) == target:
                        target_neighbors[neighbor].append(target_arrays[l][0, neighbor[0], neighbor[1]])
        for neighbor in target_neighbors:
            if any(time <= time_limit for time in target_neighbors[neighbor]):
                neighbors[(target, neighbor)] = target_neighbors[neighbor]
                
    return neighbors

def compute_neighborhoods(neighbors, target_types_at_cells, time_limit):
    neighborhoods = defaultdict(list)

    for target, neighbor in neighbors:
        time = max(neighbors[(target, neighbor)])
        if time <= time_limit:
            neighborhoods[target].append((neighbor, time))
        
    return neighborhoods

def extract_threshold_targets(neighborhoods,
                              target_types_at_cells,
                              scenarios_target_priorities,
                              scenarios_target_groups,
                              scenarios_groups_per_priority,
                              time_limit):

    scenarios_threshold_targets = [set() for scenario in scenarios_groups_per_priority]

    for target in target_types_at_cells:
        for scenario, groups_per_priority in enumerate(scenarios_groups_per_priority):
            high_priority_groups = set()
            for target_type in target_types_at_cells[target]:
                if scenarios_target_priorities[scenario][target_type] == 0:
                    high_priority_groups.add(scenarios_target_groups[scenario][target_type])
                    
            target_neighbors = deque([(target, 0.)])
            visited_neighbors = set([target])
            while target_neighbors:
                target_neighbor, target_neighbor_time = target_neighbors.pop()
                neighborhood = neighborhoods.get(target_neighbor)
                if neighborhood is not None:
                    for neighbor in neighborhood:
                        if (neighbor[0] not in visited_neighbors
                            and target_neighbor_time + neighbor[1] <= time_limit):
                            for target_type in target_types_at_cells[neighbor[0]]:
                                if (scenarios_target_priorities[scenario][target_type] == 0):
                                    high_priority_groups.add(scenarios_target_groups[scenario][target_type])
                            target_neighbors.append((neighbor[0], target_neighbor_time + neighbor[1]))
                            visited_neighbors.add(neighbor[0])
                                
            if groups_per_priority[0] == high_priority_groups:
                scenarios_threshold_targets[scenario].add(target)
                    
    return scenarios_threshold_targets

# Computation of the traverse paths

def compute_path_rank(traverse_path,
                      scenarios_visited_groups_per_priorities,
                      scenarios_path_duration,
                      max_path_length,
                      scenarios_target_priorities,
                      scenarios_target_groups,
                      scenarios_priorities,
                      scenarios_groups_per_priority,
                      high_resolution_arrays,
                      rad_longitude_array,
                      rad_latitude_array,
                      planet_radius,
                      group_weights,
                      number_weight,
                      data_weight,
                      sinuosity_weight,
                      duration_weight):
    
    ranks = [None for i in range(len(scenarios_visited_groups_per_priorities))]
    if scenarios_visited_groups_per_priorities == ranks:
        return ranks
    
    path_length = 0.
    previous_longitude = rad_longitude_array[traverse_path[0][0][0]]
    previous_latitude = rad_latitude_array[traverse_path[0][0][0]]
    current_longitude = rad_longitude_array[traverse_path[-1][0][0]]
    current_latitude = rad_latitude_array[traverse_path[-1][0][0]]
    path_straight_length = haversine_distance_math(current_longitude, current_latitude,
                                                   previous_longitude, previous_latitude, planet_radius)
    number_highres_data = 0

    for target in traverse_path:
        current_longitude = rad_longitude_array[target[0][0]]
        current_latitude = rad_latitude_array[target[0][0]]
        path_length += haversine_distance_math(current_longitude, current_latitude,
                                               previous_longitude, previous_latitude, planet_radius)
        previous_longitude = current_longitude
        previous_latitude = current_latitude

        data_ratio = 0.
        for array in high_resolution_arrays:
            if array[target[0][0]] == 1:
                data_ratio += 1
        if len(high_resolution_arrays) > 0:
            number_highres_data += data_ratio/len(high_resolution_arrays)
    
    data_rank = 0.
    if len(traverse_path) > 0:
        data_rank = number_highres_data/len(traverse_path)
    path_sinuosity = 0.
    if path_length > 0.:
        path_sinuosity = path_straight_length/path_length
    
    for scenario in range(len(scenarios_visited_groups_per_priorities)):
        if scenarios_visited_groups_per_priorities[scenario] is not None:
            duration_rank = 0.
            if max_path_length > 0.:
                duration_rank = 1 - scenarios_path_duration[scenario]/max_path_length
            ranks[scenario] = data_weight*data_rank + sinuosity_weight*path_sinuosity + duration_weight*duration_rank
            for path_groups, target_groups, group_weight in zip(scenarios_visited_groups_per_priorities[scenario],
                                                                scenarios_groups_per_priority[scenario],
                                                                group_weights[scenario]):
                group_rank = 0.
                if len(target_groups) > 0:
                    group_rank = len(set(path_groups))/len(target_groups)
                number_rank = 0.
                if len(path_groups) > 0:
                    number_rank = 1 - len(set(path_groups))/len(path_groups)
                ranks[scenario] += group_weight*group_rank + number_weight*number_rank

    return ranks

def are_all_high_priority_in_path(traverse_path,
                                  scenarios_groups_per_priority,
                                  scenarios_target_priorities,
                                  scenarios_target_groups):
    scenarios_visited_groups_per_priorities = [None for scenario in scenarios_groups_per_priority]
    scenarios_path_duration = [None for scenario in scenarios_groups_per_priority]
    for scenario in range(len(scenarios_visited_groups_per_priorities)):
        temp = [list() for groups in scenarios_groups_per_priority[scenario]]
        path_duration = 0
        for target in traverse_path:
            if scenarios_target_priorities[scenario][target[0][1]] is not None:
                if set(temp[0]) != scenarios_groups_per_priority[scenario][0]:
                    path_duration = target[1]
                temp[scenarios_target_priorities[scenario][target[0][1]]].append(scenarios_target_groups[scenario][target[0][1]])
        if set(temp[0]) == scenarios_groups_per_priority[scenario][0]:
            scenarios_visited_groups_per_priorities[scenario] = temp
            scenarios_path_duration[scenario] = path_duration
            
    return scenarios_visited_groups_per_priorities, scenarios_path_duration

def check_path_validity(traverse_path,
                        new_target,
                        max_path_duration):
    new_path_duration = traverse_path[-1][1] + new_target[1]
    is_new_target_in_path = False
    i = 0
    while i < len(traverse_path) and is_new_target_in_path == False:
        if traverse_path[i][0][0] == new_target[0]:
            is_new_target_in_path = True
        i += 1
    if (new_path_duration <= max_path_duration and
        is_new_target_in_path == False):
        return True
    return False

def compute_traverse_paths(threshold_targets,
                           neighborhoods,
                           target_types_at_cells,
                           max_path_length,
                           scenarios_target_priorities,
                           scenarios_target_groups,
                           scenarios_priorities,
                           scenarios_groups_per_priority,
                           high_resolution_arrays,
                           rad_longitude_array,
                           rad_latitude_array,
                           planet_radius,
                           group_weights,
                           number_weight,
                           data_weight,
                           sinuosity_weight,
                           duration_weight):
    
    traverse_paths = defaultdict(list)
    
    for target in threshold_targets:

        initial_path = [((target, target_type), 0.) for target_type in target_types_at_cells[target]]
        open_paths = deque([initial_path])
        final_paths = [[] for scenario in range(len(scenarios_priorities))]

        while open_paths:

            current_path = open_paths.pop()
            last_cell = current_path[-1][0][0]
            neighborhood = neighborhoods.get(last_cell)

            was_neighbor_added = False
            if neighborhood is not None:
                for neighbor in neighborhood:
                    if check_path_validity(current_path,
                                           neighbor,
                                           max_path_length) == True:
                        new_path = list(current_path)
                        for target_type in target_types_at_cells[neighbor[0]]:
                            new_path.append(((neighbor[0], target_type),
                                              current_path[-1][1] + neighbor[1]))
                        open_paths.append(new_path)
                        was_neighbor_added = True
            if was_neighbor_added == False:
                scenarios_visited_groups_per_priorities, scenarios_path_duration = are_all_high_priority_in_path(current_path,
                                                                                                                 scenarios_groups_per_priority,
                                                                                                                 scenarios_target_priorities,
                                                                                                                 scenarios_target_groups)

                path_ranks = compute_path_rank(current_path,
                                               scenarios_visited_groups_per_priorities,
                                               scenarios_path_duration,
                                               max_path_length,
                                               scenarios_target_priorities,
                                               scenarios_target_groups,
                                               scenarios_priorities,
                                               scenarios_groups_per_priority,
                                               high_resolution_arrays,
                                               rad_longitude_array,
                                               rad_latitude_array,
                                               planet_radius,
                                               group_weights,
                                               number_weight,
                                               data_weight,
                                               sinuosity_weight,
                                               duration_weight)

                for scenario, rank in enumerate(path_ranks):
                    if rank is not None:
                        final_paths[scenario].append((rank, scenarios_path_duration[scenario], current_path))
                    
        for scenario, paths in enumerate(final_paths):
            if paths:
                sorted_paths = sorted(paths)
                traverse_paths[tuple(sorted_paths[-1][2])].append((sorted_paths[-1][0],
                                                                   sorted_paths[-1][1],
                                                                   scenario))
                        
    return traverse_paths

def save_paths_to_csv_file(file_path, paths_dict):
    with open(file_path, "w") as output:
        writer = csv.writer(output)
        writer.writerows(paths_dict.items())

def read_paths_from_csv_file(file_path):
    traverse_paths = dict()

    with open(file_path, "r") as output:
        reader = csv.reader(output)
        for row in reader:
            path = make_tuple(row[0])
            traverse_paths[path] = make_tuple(row[1])

    return traverse_paths