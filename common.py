#!/usr/bin/env python

"""Common functions."""

import numpy as np

detector_length = 90.0
detector_width = 47.5
detector_height = 40.0
xx = (0, detector_width)
yy = (-detector_height/2.0, detector_height/2.0)
zz = (0, detector_length)
x_ = (3.0, detector_width-3.0)
y_ = (-detector_height/2.0+4.0, detector_height/2.0-4.0)
z_ = (6.0, detector_length-4.0)

def inside_tpc(x, y, z):
    """Returns True if x, y, and z are inside active TPC volume."""
    x = round(x, 9)
    y = round(y, 9)
    z = round(z, 9)
    if (x >= xx[0] and x <= xx[1] and y >= yy[0] and y <= yy[1] and
        z >= zz[0] and z <= zz[1]):
        return True
    else:
        return False

def inside_fiducial(x, y, z):
    """Returns True if x, y, and z are inside fiducial volume."""
    x = round(x, 9)
    y = round(y, 9)
    z = round(z, 9)
    if (x >= x_[0] and x <= x_[1] and y >= y_[0] and y <= y_[1] and
        z >= z_[0] and z <= z_[1]):
        return True
    else:
        return False

def tpc_exit_point(p1, p2):
    """
    Returns the point at which the shower exits the active TPC volume.

    p1 is the end point of the charged pion (where the pion charge-
    exchange occurs).
    p2 is the end point of the photon from the neutral pion decay
    (where the photon pair converts).

    """
    e0 = (xx[0], yy[0], zz[0])
    e1 = (xx[1], yy[1], zz[1])
    p1 = np.array(p1, dtype=np.float64)
    p2 = np.array(p2, dtype=np.float64)
    v = p2 - p1
    v /= np.sqrt(np.dot(v, v))
    parameters = np.array([], dtype=np.float64)
    for i in xrange(3):
        if v[i] == 0:
            continue
        parameters = np.append(parameters, (e0[i] - p1[i]) / v[i])
        parameters = np.append(parameters, (e1[i] - p1[i]) / v[i])
    #parameters = np.append(parameters, (xx[0] - p1[0]) / v[0])
    #parameters = np.append(parameters, (xx[1] - p1[0]) / v[0])
    #parameters = np.append(parameters, (yy[0] - p1[1]) / v[1])
    #parameters = np.append(parameters, (yy[1] - p1[1]) / v[1])
    #parameters = np.append(parameters, (zz[0] - p1[2]) / v[2])
    #parameters = np.append(parameters, (zz[1] - p1[2]) / v[2])
    parameters = parameters[parameters > 0]
    tpc_edge = np.array([ 0, 0, 0 ], dtype=np.float64)
    for i in xrange(parameters.size):
        x = p1[0] + parameters[i]*v[0]
        y = p1[1] + parameters[i]*v[1]
        z = p1[2] + parameters[i]*v[2]
        if inside_tpc(x, y, z):
            tpc_edge = np.array([ x, y, z ], dtype=np.float64)
    return tpc_edge

def distance_from_tpc_edge(p1, p2):
    """
    Returns the distance from photon conversion point to TPC edge along
    the direction of photon propagation.

    p1 is the end point of the charged pion (where the pion charge-
    exchange occurs).
    p2 is the end point of the photon from the neutral pion decay
    (where the photon pair converts).

    """
    p1 = np.array(p1, dtype=np.float64)
    p2 = np.array(p2, dtype=np.float64)
    tpc_edge = tpc_exit_point(p1, p2)
    tpc_edge_vector = tpc_edge - p2
    distance = np.sqrt(np.dot(tpc_edge_vector, tpc_edge_vector))
    return distance

def distance_from_tpc_z_axis(p):
    """
    Returns the distance from point p z-axis of the TPC.

    """
    p = np.array(p, dtype=np.float64)
    distance = np.sqrt(np.square(p[0]-detector_width/2.0) + np.square(p[1]))
    return distance

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

def cosine_between(v1, v2):
    """ Returns the cosine between vectors 'v1' and 'v2'::

            >>> cosine_between((1, 0, 0), (0, 1, 0))
            0.0
            >>> cosine_between((1, 0, 0), (1, 0, 0))
            1.0
            >>> cosine_between((1, 0, 0), (-1, 0, 0))
            -1.0
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    cosine = np.dot(v1_u, v2_u)
    return cosine

def pion_mass(e1, e2, cos):
    """ Returns the mass of pi0 given gamma energies and angle """
    return np.sqrt(2*e1*e2*(1-cos))

