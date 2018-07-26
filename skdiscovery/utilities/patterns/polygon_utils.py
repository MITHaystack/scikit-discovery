import numpy as np
from shapely.geometry import Polygon, Point
from collections import OrderedDict

def shoelaceArea(in_vertices):
    x = in_vertices[:,0]
    y = in_vertices[:,1]
    return 0.5 *(np.sum(x * np.roll(y,shift=-1)) - np.sum(np.roll(x,shift=-1) * y))


def parseBasemapShape(aquifers, aquifers_info):
    polygon_data = []
    test_list = []

    for index,(aquifer,info) in enumerate(zip(aquifers,aquifers_info)):
        if shoelaceArea(np.array(aquifer)) < 0:
            new_data = OrderedDict()
            new_data['shell'] = aquifer
            new_data['info'] = info
            new_data['holes'] = []

            polygon_data.append(new_data)

        else:
            polygon_data[-1]['holes'].append(aquifer)


    for data in polygon_data:
        data['polygon'] = Polygon(shell=data['shell'],holes=data['holes'])

    return polygon_data



def nearestEdgeDistance(x,y,poly):

    point = Point(x,y)

    ext_dist = poly.exterior.distance(point)

    if len(poly.interiors) > 0:
        int_dist = np.min([interior.distance(point) for interior in poly.interiors])

        return np.min([ext_dist, int_dist])
    else:
        return ext_dist


def findPolygon(in_data, in_point):
    result_num = None
    for index, data in enumerate(in_data):
        if data['polygon'].contains(in_point):
            if result_num == None:
                result_num = index
            else:
                raise RuntimeError("Multiple polygons contains point")

    if result_num == None:
        return -1
    return result_num

def getInfo(row, key, fill, polygon_data):
    try:
        return polygon_data[int(row['ShapeIndex'])]['info'][key]
    except KeyError:
        return fill


def findClosestPolygonDistance(x,y,polygon_data):
    min_dist = np.inf
    shape_index = -1
    point = Point(x,y)
    for index, data in enumerate(polygon_data):
        if not data['polygon'].contains(point) and data['info']['AQ_CODE'] != 999:
            new_distance = data['polygon'].distance(point)
            if  new_distance < min_dist:
                min_dist = new_distance
                shape_index = index

    return min_dist, shape_index
