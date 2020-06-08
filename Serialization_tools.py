import numpy as np
import Arithmetic as arithm

def deconstruct(points):
    tmp_points = []
    for P in points:
        tmp_x = [coord.x for coord in P[0].x]
        tmp_y = [coord.x for coord in P[1].x]
        tmp_points.append([tmp_x,tmp_y])
    return np.array(tmp_points)

def reconstruct(points,p):
    tmp_points = [np.array([arithm.ZPIfield(P[0][0], P[0][1], p), arithm.ZPIfield(P[1][0],P[1][1], p)]) for P in points]
    return np.array(tmp_points)

def serialize_points(path,points):
    to_save = deconstruct(points)
    with open(path, 'wb+') as f:
        np.save(f, to_save)

def deserialize_points(path,p):
    with open(path,'rb') as f:
        points = np.load(f)
    return reconstruct(points,p)