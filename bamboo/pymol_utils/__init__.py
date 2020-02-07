import os
import numpy

from itertools import cycle

from bamboo.pymol_utils.shapes import *


class PymolColourPalettes:


    set1 = ['carbon','cyan','lightmagenta','yellow','salmon','hydrogen','slate','orange']
    set2 = ['lime','deepteal','hotpink','yelloworange','violetpurple','grey70','marine','olive']
    set3 = ['smudge','teal','dirtyviolet','wheat','deepsalmon','lightpink','aquamarine','paleyellow']
    set4 = ['limegreen','skyblue','warmpink','limon','violet','bluewhite','greencyan','sand']
    set5 = ['forest','lightteal','darksalmon','splitpea','raspberry','grey50','deepblue','brown']


class _PymolScript(object):


    _colours = PymolColourPalettes.set1

    _styles = ['cartoon', 'sticks', 'ribbon', 'dots', 'spheres', 'surface', 'mesh']

    def __init__(self, pretty_models=True, pretty_maps=True):

        # Previous object manipulated
        self.prev = None

        self.lines = []
        self._add(self._cmd_import_pml)
        self._add(self._cmd_import_cgo)

        if pretty_models: self.pretty_models()
        if pretty_maps:   self.pretty_maps()

        self.colour_palette(self._colours)

    # ----------------------------------------------------------------------- #
    # Access/Modify Object Properties
    # ----------------------------------------------------------------------- #

    def colour_palette(self, colours):
        self._colours = colours
        self._colour_cycle = cycle(self._colours)

    def get_colour_palette(self):
        return self._colours

    def get_colour(self):
        return self._colour_cycle.next()

    # ----------------------------------------------------------------------- #
    # Append to/write script
    # ----------------------------------------------------------------------- #

    def _add(self, thing):
        if isinstance(thing, list):     [self._add(t) for t in thing]
        elif isinstance(thing, str):    self.lines.append(thing)
        else: raise Exception('Invalid type: {}'.format(type(thing)))

    def write_script(self, f_name, overwrite=False):
        if (overwrite is True) and os.path.exists(f_name):
            os.remove(f_name)
        assert not os.path.exists(f_name)
        assert f_name.endswith(self._file_type)
        with open(f_name, 'w') as fh:
            fh.write('\n'.join(self.lines)+'\n')

    # ----------------------------------------------------------------------- #
    # Miscellanous Commands
    # ----------------------------------------------------------------------- #

    def ray_trace(self, width=1024, height=786):
        self._add(self._cmd_ray.format(width=width, height=height))

    def quit(self):
        self._add(self._cmd_quit)

    # ----------------------------------------------------------------------- #
    # Modify the settings in PyMOL session
    # ----------------------------------------------------------------------- #

    def set(self, setting, *args, **kwargs):
        self._add(self._cmd_set.format(setting=setting,
                                       args=', '.join([repr(a) for a in args]),
                                       kwargs=', '.join(['{}={}'.format(k, repr(kwargs[k])) for k in kwargs.keys()])))

    def set_normalise_maps(self, value=True):
        self.set("normalize_ccp4_maps", int(value))

    # ----------------------------------------------------------------------- #
    # Load and fetch structures and maps
    # ----------------------------------------------------------------------- #

    def fetch_pdb(self, pdb_code):
        pass

    def fetch_map(self, pdb_code, map_type='2fofc'):
        pass

    def load_pdb(self, f_name, obj=None):
        """Load pdb file f_name"""
        if obj is None: obj = os.path.basename(os.path.splitext(f_name)[0])
        self._add(self._cmd_load_basic.format(f_name=f_name, obj=obj))
        return obj

    def load_map(self, f_name, obj=None):
        if obj is None: obj = os.path.basename(os.path.splitext(f_name)[0])
        self._add(self._cmd_load_basic.format(f_name=f_name, obj=obj))
        return obj

    def make_mesh(self, obj, mesh_suffix='.mesh', contour_level=1, colour=None):
        mesh_name = obj+mesh_suffix
        self._add(self._cmd_isomesh.format(obj=mesh_name, map=obj, contour_level=contour_level))
        if colour is not None: self.colour(obj=mesh_name, colour=colour)
        return mesh_name

    # ----------------------------------------------------------------------- #
    # Quick addition of Shape objects
    # ----------------------------------------------------------------------- #

    def add_shape(self, shape, obj):
        self._add(shape.as_cmd(name=obj))

    # ----------------------------------------------------------------------- #
    # Multi-command defaults
    # ----------------------------------------------------------------------- #

    def pretty_models(self):
        self.set("cartoon_side_chain_helper", 1)

    def pretty_maps(self):
        pass

    # ----------------------------------------------------------------------- #
    # Visualisation
    # ----------------------------------------------------------------------- #

    def repr_as(self, obj, style='cartoon'):
        assert style in self._styles
        self._add(self._cmd_show_as.format(style=style, obj=obj))
        return style

    def repr_show(self, obj, style='cartoon'):
        assert style in self._styles
        self._add(self._cmd_show.format(style=style, obj=obj))
        return style

    def repr_hide(self, obj, style='everything'):
        self._add(self._cmd_hide.format(style=style, obj=obj))
        return style

    def colour(self, obj, colour=None):
        if colour is None: colour = self.get_colour()
        self._add(self._cmd_colour.format(colour=colour, obj=obj))
        return colour

    # ----------------------------------------------------------------------- #

class PythonScript(_PymolScript):


    _file_type = '.py'

    _cmd_import_pml = 'from pymol import cmd'
    _cmd_import_cgo = 'from pymol.cgo import *'

    _cmd_quit       = 'cmd.quit()'

    _cmd_ray        = 'cmd.ray({width},{height})'

    _cmd_load_basic = 'cmd.load("{f_name}","{obj}")'
    _cmd_isomesh    = 'cmd.isomesh("{obj}", "{map}", {contour_level})'

    _cmd_show_as    = 'cmd.show_as("{style}", "{obj}")'
    _cmd_show       = 'cmd.show("{style}", "{obj}")'
    _cmd_hide       = 'cmd.hide("{style}", "{obj}")'
    _cmd_colour     = 'cmd.color("{colour}", "{obj}")'
    _cmd_set        = 'cmd.set("{setting}", {args}, {kwargs})'


