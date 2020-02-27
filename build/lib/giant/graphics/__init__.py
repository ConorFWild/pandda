import math

# Code taken from Dimple (Marcin)

def calculate_view_quaternion(p1, p2):

    if p1 is None or p2 is None:
        return (0., 0., 0., 1.)
    d = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
    length = math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
    d = (d[0]/length, d[1]/length, d[2]/length)
    # ref vector: 0 0 -1 (?)
    # cross product (a2 b3 - a3 b2, ..., ...)
    prod = (d[1], -d[0], 0)
    quat = (prod[0], prod[1], prod[2], 1-d[2])
    qlen = math.sqrt(sum(a*a for a in quat))
    return (quat[0]/qlen, quat[1]/qlen, quat[2]/qlen, quat[3]/qlen)

def multiply_quaternions(q1, q2):
    x, y, z, w = q1
    ax, ay, az, aw = q2
    return (w*ax + x*aw + y*az - z*ay,
            w*ay + y*aw + z*ax - x*az,
            w*az + z*aw + x*ay - y*ax,
            w*aw - x*ax - y*ay - z*az)
