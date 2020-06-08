import numpy as np

#Function returns True if given array contains given element.
#Otherwise returns False
def isInArray(array, x):
    for item in array:
        if item[0]==x[0] and item[1]==x[1]:
            return True
    return False

#Returns set described in Velu algorithm from given subgroup of elliptic curve
def Velu_set(points,el):
    rank_two = np.array([P for P in points if el.m_mltpl(2,P)==0])
    G2set = []
    tmp_R = []
    R_plus = []
    R_minus = []
    for P in points:
        if not isInArray(rank_two, P):
            tmp_R.append(P)
        else:
            G2set.append(P)
    for P in tmp_R:
        minus_P = el.opposite(P)
        if not (isInArray(R_plus, P) or isInArray(R_minus, P)):
            R_plus.append(P)
            R_minus.append(minus_P)
    if len(G2set)>0:
        return np.concatenate((np.array(R_plus), np.array(G2set)), axis=0)
    return np.array(R_plus)

#Returns coefficients of image curve and other values desribed in Velu algorithm
def Velu_values(V,el):
    v = el.const(0)
    w = el.const(0)
    aux_values = []
    for P in V:
        gqx = el.const(3) * P[0] * P[0] + el.const(2) * el.coeffs[1] * P[0] + el.coeffs[3] - el.coeffs[0] * P[1]
        gqy = el.const(0) - el.const(2) * P[1] - el.coeffs[0] * P[0] - el.coeffs[2]
        if el.m_mltpl(2,P) == 0:
            vq = gqx
        else:
            vq = el.const(2) * gqx - el.coeffs[0] * gqy
        uq = gqy * gqy
        aux_values.append({'gqx': gqx, 'gqy': gqy, 'vq': vq, 'uq': uq, 'xq': P[0], 'yq': P[1]})
        v = v + vq
        w = w + uq + P[0] * vq
    new_coeffs = [el.coeffs[0], el.coeffs[1], el.coeffs[2], el.coeffs[3] - el.const(5) * v,
                  el.coeffs[4] - (el.coeffs[0] * el.coeffs[0] + el.const(4) * el.coeffs[1]) * v - el.const(7) * w]
    return np.array(new_coeffs),np.array(aux_values)

#returns x-coordinate rational function of Velu isogeny
def Velu_x(aux_values,el):
    return lambda x,y : x + np.sum(np.array([val['vq']*((x-val['xq']).inv()) + val['uq']*(((x-val['xq'])*\
                        (x-val['xq'])).inv()) for val in aux_values]))

#returns y-coordinate rational function of Velu isogeny
def Velu_y(aux_values,el):
    return lambda x,y : y - np.sum(np.array([val['uq']*(el.const(2)*y+el.coeffs[0]*x+el.coeffs[2])*\
                        (((x-val['xq'])*(x-val['xq'])*(x-val['xq'])).inv()) + val['vq']*(el.coeffs[0]*(x-val['xq'])+y-val['yq'])*\
                        (((x-val['xq'])*(x-val['xq'])).inv()) + (el.coeffs[0]*val['uq']-val['gqx']*val['gqy'])*\
                        (((x-val['xq'])*(x-val['xq'])).inv()) for val in aux_values]))