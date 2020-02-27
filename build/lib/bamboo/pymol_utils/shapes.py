try:
    from pymol import cmd
    from pymol.cgo import *
except:
    pass

class Shape(object):
    def draw(self, name):
        print 'Loading %s to %s ' % (self.SHAPE, name)
        exec(self.as_cmd(name))
    def as_cmd(self, name):
        return 'cmd.load_cgo(%s, "%s")' % (self.as_string(), name)
    def as_string(self):
        raise Exception('Not Defined')

class Cylinder(Shape):
    SHAPE = 'Cylinder'
    def __init__(self, start, end, colors=None, radius=None):
        if (colors is None):    colors = [[1,0,0]]*2
        if (radius is None):     radius  = 0.1
        self.start = start
        self.end = end
        self.colors = colors
        self.radius = radius

    def as_string(self):
        s = "["
        s += "CYLINDER," + " %.6g, %.6g, %.6g," % tuple(self.start)
        s += " %.6g, %.6g, %.6g," % tuple(self.end)
        s += " %.6g," % self.radius
        s += " %.6g, %.6g, %.6g," % tuple(self.colors[0])
        s += " %.6g, %.6g, %.6g" % tuple(self.colors[1])
        s += "]"
        return s

class Sphere(Shape):
    SHAPE = 'Sphere'
    def __init__(self, centre, radius=None):
        if (radius is None):    radius  = 1
        self.centre = centre
        self.radius = radius
#        self.color = color
    def as_string(self):
        s = "["
        s += "SPHERE,"
        s += " %.6g, %.6g, %.6g," % tuple(self.centre)
        s += " %.6g" % self.radius
        s += "]"
        return s

