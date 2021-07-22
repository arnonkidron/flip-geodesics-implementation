from math import acos, cos, sin, sqrt, pi
import numpy as np


################################################
# the cosine theorem for angles and side lengths
################################################
def get_angle(len1, len2, len3):
    """
    Input: the lengths of 3 sides of a triangle
    Output: the angle between len1, len2, opposing len3, in radians
    """
    return acos((len1 * len1 + len2 * len2 - len3 * len3) / (2.0 * len1 * len2))


def get_side_length(a, b, angle):
    """
    Input: the lengths of 2 sides
           the angle between them
    Output: the length of the other side

    """
    return sqrt(a * a + b * b - 2 * a * b * cos(angle))


def get_angle_between(vec1, vec2):
    return acos(np.dot(vec1, vec2))


def is_reflex(angle):
    return angle >= pi - 1e-09


#################################
# general rotation around an axis
#################################

def rotate(v, theta, n):
    """
    :arg v: a vector
    :arg theta: an angle
    :arg n: the axis to rotate around
    :return: the vector obtained by rotating v around n by theta
    """
    v = np.array(v)
    n = np.array(n)
    n = n / np.linalg.norm(n)
    cost = cos(theta)
    sint = sin(theta)
    return v * cost \
        + n * np.dot(n, v) * (1 - cost) \
        + np.cross(n, v) * sint


def turn(v, theta, towards):
    """
    :arg v: a vector
    :arg theta: an angle
    :arg towards: another vector
    :return: rotate v by theta towards the other vector
    """
    v = np.array(v)
    towards = np.array(towards)
    # result = rotate(v, theta, np.cross(v, towards))
    # return result
    return rotate(v, theta, np.cross(v, towards))


###########################

def get_closest_point(start_A, vec_A, start_B, vec_B):
    """
    Input: 2 lines A,B, determined by a start point and a unit vector
    Output: the point on line A that is closest to line B
    source: https://math.stackexchange.com/questions/1993953/closest-points-between-two-lines
    """
    vec_C = np.cross(vec_A, vec_B)
    rhs = start_B - start_A
    lhs = np.array([vec_A, vec_B, vec_C]).T
    t1, t2, t3 = np.linalg.solve(lhs, rhs)
    return start_A + vec_A * t1
