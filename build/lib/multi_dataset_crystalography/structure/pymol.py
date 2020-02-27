import os

from iotbx.pdb.hierarchy import input as ih

from bamboo.pymol_utils import PythonScript
from giant.structure.select import get_select_function
from giant.structure.formatting import PymolSelection, ShortLabeller

def auto_chain_images(structure_filename, output_prefix, selection='protein', style='cartoon', **kw_args):
    filter_func = get_select_function(selection)
    return selection_images(structure_filename  = structure_filename,
                            output_prefix       = output_prefix+'-chain_',
                            selections  = [PymolSelection.format(c) for c in filter_func(ih(structure_filename).hierarchy).chains()],
                            labels      = [ShortLabeller.format(c)  for c in filter_func(ih(structure_filename).hierarchy).chains()],
                            style = style, **kw_args)

def auto_residue_images(structure_filename, output_prefix, selection='protein', style='sticks', **kw_args):
    filter_func = get_select_function(selection)
    return selection_images(structure_filename  = structure_filename,
                            output_prefix       = output_prefix+'-residue_',
                            selections  = [PymolSelection.format(r) for r in filter_func(ih(structure_filename).hierarchy).residue_groups()],
                            labels      = [ShortLabeller.format(r)  for r in filter_func(ih(structure_filename).hierarchy).residue_groups()],
                            style = style, **kw_args)

def selection_images(structure_filename,
                     output_prefix,
                     selections,
                     labels = None,
                     style  = 'sticks',
                     colour_by = None,
                     hide_rest = True,
                     ray_trace = True,
                     settings = [],
                     run_script = True,
                     delete_script  = True,
                     width  = 800,
                     height = 600,
                    ):

    if labels is None:
        labels = map(str,range(1, len(selections)+1))

    # Create script object
    s = PythonScript(pretty_models=False, pretty_maps=False)
    s.custom('bg_color', 'white')
    s.set('opaque_background', 0)
    s.set('ray_opaque_background', 0)
    s.set('orthoscopic', 1)
    # Apply custom global commands
    for cmd in settings:
        s.set(*cmd)
    # Read in structure and hide all atoms
    s.load_pdb(f_name=structure_filename, obj='input_structure')
    # Set styles
    style = style.split('+')

    png_filenames = []

    for i, selection in enumerate(selections):
        # Reset the structure to all look the same
        s.show_as(obj='all', style=style[0])
        s.colour(obj='all', colour="grey90")
        # Apply custom views to the selection
        s.select(obj='sele', selection=selection)
        for sty in style:
            s.show(obj='sele', style=sty)
        s.orient(obj='sele')
        s.zoom(obj='sele', buffer=0.0, complete=1)
        if colour_by == 'bfactor':
            s.custom('spectrum', expression='b', selection='sele and (b>0)')
        else:
            s.colour_by_element(obj='sele', carbon_colour='green')
        if hide_rest:
            s.hide(obj='not sele', style='everything')
        if ray_trace:
            s.ray(height=height, width=width)
        png_name = s.png(f_name=output_prefix+labels[i]+'.png')
        png_filenames.append(png_name)

    f_name = s.write_script(output_prefix+'.py')
    l_name = f_name.replace('.py', '.log')
    assert not os.path.exists(l_name)

    if run_script is True:
        s.run(f_name)

    if delete_script is True:
        if os.path.exists(f_name):
            os.remove(f_name)
        if os.path.exists(l_name):
            os.remove(l_name)

    return png_filenames
