import numpy as np
import Arithmetic as ar

def deconstruct(points):
    tmp_points = []
    for P in points:
        tmp_x = []
        tmp_y = []
        for coord in P[0].x:
            tmp_x.append(coord.x)
        for coord in P[1].x:
            tmp_y.append(coord.x)
        tmp_points.append([tmp_x,tmp_y])
    return np.array(tmp_points)

def reconstruct(points,p):
    tmp_points = []
    for P in points:
        tmp_points.append(np.array([ar.ZPIfield(P[0][0], P[0][1], p), ar.ZPIfield(P[1][0],P[1][1], p)]))
    return np.array(tmp_points)

def serialize_points(path,points):
    to_save = deconstruct(points)
    with open(path, 'wb+') as f:
        np.save(f, to_save)

def deserialize_points(path,p):
    with open(path,'rb') as f:
        points = np.load(f)
    return reconstruct(points,p)