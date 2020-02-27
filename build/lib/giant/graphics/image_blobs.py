
import os
from pandda.html import PANDDA_HTML_ENV
from giant.graphics import calculate_view_quaternion

# Get the template to be filled in
template = PANDDA_HTML_ENV.get_template('ccp4mg-pic.py')

point = [13.32,36.95,26.02]
com = [16.16,37.38,18.6]
orientation = calculate_view_quaternion(p1=com, p2=point)

ccp4mg_script = template.render({
                                    'view'  :{    'camera_centre' : [-1.0*p for p in point],
                                                  'orientation'   : orientation
                                                },
                                    'mol'   :{    'file'  : 'aligned_structure',
                                                  'path'  : './aligned.pdb',
                                                  'name'  : 'aligned_structure'
                                                },
                                    'map'   :{    'file'    : 'sampled_map',
                                                  'path'    : './sampled.ccp4',
                                                  'name'    : 'sampled_map',
                                                  'contour' : 1
                                                },
                                    'diff_map' :{ 'file'    : 'diff_map',
                                                  'path'    : './z_map.ccp4',
                                                  'name'    : 'diff_map',
                                                  'pos-contour' : 3,
                                                  'neg-contour' : -3
                                                }
                                })

print '################################'
print ccp4mg_script
print '################################'

outfile = './ccpm4_script.py'
assert not os.path.exists(outfile)

open(outfile, 'w').write(ccp4mg_script)











