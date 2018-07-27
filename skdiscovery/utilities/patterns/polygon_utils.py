import numpy as np
from shapely.geometry import Polygon, Point
from collections import OrderedDict

def shoelaceArea(in_vertices):
    """
    Determine the area of a polygon using the shoelace method

    https://en.wikipedia.org/wiki/Shoelace_formula

    @param in_vertices: The vertices of a polygon. 2d Array where the first column is the
                        x coordinates and the second column is the y coordinates

    @return: Area of the polygon
    """
    x = in_vertices[:,0]
    y = in_vertices[:,1]
    return 0.5 *(np.sum(x * np.roll(y,shift=-1)) - np.sum(np.roll(x,shift=-1) * y))


def parseBasemapShape(aquifers, aquifers_info):
    """
    Create shapely polygons from shapefile read in with basemap

    @param aquifers: Data read in shapefile from basemap
    @param aquifers_info: Metadata read from shapefile from basemap

    @return: Dictionary containing information about shapes and shapely polygon of shapefile data
    """
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
    """
    Determine the distance to the closest edge of a polygon

    @param x: x coordinate
    @param y: y coordinate
    @param poly: Shapely polygon

    @return distance from x,y to nearest edge of the polygon
    """
    point = Point(x,y)

    ext_dist = poly.exterior.distance(point)

    if len(poly.interiors) > 0:
        int_dist = np.min([interior.distance(point) for interior in poly.interiors])

        return np.min([ext_dist, int_dist])
    else:
        return ext_dist


def findPolygon(in_data, in_point):
    """
    Find the polygon that a point resides in

    @param in_data: Input data containing polygons as read in by parseBasemapShape
    @param in_point: Shapely point

    @return: Index of shape in in_data that contains in_point
    """
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
    """
    Retrieve information from polygon data:

    @param row: Container with key 'ShapeIndex'
    @param key: Key of data to retrieve from polygon_data element
    @param fill: Value to return if key does not exist in polygon_data element
    @param polygon_data: Polygon data as read in by parseBasemapShape
    """
    try:
        return polygon_data[int(row['ShapeIndex'])]['info'][key]
    except KeyError:
        return fill


def findClosestPolygonDistance(x,y,polygon_data):
    """
    Find the distance to the closest polygon

    @param x: x coordinate
    @param y: y coordinate
    @param polygon_data: Polygon data as read in by parseBasemapShape
    @return Distance from x, y to the closest polygon polygon_data
    """
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
