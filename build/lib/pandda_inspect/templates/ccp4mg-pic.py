
#
# Open with the command
# ccp4mg -norestore -picture testpic.mgpic.py
#
# Save image and quit with
# ccp4mg -norestore -picture testpic.mgpic.py -R file-name.png -RO '{"size":"1000x1000"}' -quit
#

# SET THE VIEW
view = View(
    centre_xyz = {{ view['camera_centre'] }},
    orientation = {{ view['orientation'] }},
#    scale = 'auto',
    scale = 30,
    fog_enabled = True,
    fog_near = -0.0,
    fog_far = 20.0,
    slab_enabled = 1,
    slab_width = 10,
    slab_offset = 0,
    clip_cap = False,
    customClipPlanes = [],
)

{% if mol %}
# OPEN MOL
MolData(
    filename =  [   'FULLPATH',
                    '{{ mol['path'] }}',
                    '{{ mol['path'] }}'
                ],
    name = '{{ mol['name'] }}',
)

# DISPLAY MOL
MolDisp(
    selection_parameters = {
                            'select' : 'cid',
                            'cid' : 'all',
                           },
    colour_parameters = {
                            'colour_mode' : 'atomtype'
                        },
)
{% endif %}

{% if map %}
# OPEN MAP
MapData(
    filename = [    'FULLPATH',
                    '{{ map['path'] }}',
                    '{{ map['path'] }}'
                ],
    name = '{{ map['name'] }}',
    expanded = 1,
    filetype = 'MAP',
)

{% for contour in map['contour'] %}
# DISPLAY DIFF MAP
MapDisp(
    opacity = 1.0,
    style = 'surface_style_solid_chickenwire',
    radius = 7,
#    clip_mode = 'ON',
    surface_style = 'surface_solid_chickenwire',
    contour_level = {{ contour }},
    contour_scale = 'absolute',
    colour = 'blue',
)
{% endfor %}
{% endif %}

{% if diff_map %}
# OPEN MAP
MapData(
    filename = [    'FULLPATH',
                    '{{ diff_map['path'] }}',
                    '{{ diff_map['path'] }}'
                ],
    name = '{{ diff_map['name'] }}',
    expanded = 1,
    filetype = 'MAP',
)

{% if diff_map['pos-contour'] %}
{% for contour in diff_map['pos-contour'] %}
# DISPLAY DIFF MAP
MapDisp(
    opacity = 1.0,
    style = 'surface_style_solid_chickenwire',
    radius = 7,
#    clip_mode = 'ON',
    surface_style = 'surface_solid_chickenwire',
    contour_level = {{ contour }},
    contour_scale = 'absolute',
    colour = 'green',
)
{% endfor %}
{% endif %}

{% if diff_map['neg-contour'] %}
{% for contour in diff_map['neg-contour'] %}
MapDisp(
    opacity = 1.0,
    style = 'surface_style_solid_chickenwire',
    radius = 7,
#    clip_mode = 'ON',
    surface_style = 'surface_solid_chickenwire',
    contour_level = {{ diff_map['neg-contour'] }},
    contour_scale = 'absolute',
    colour = 'red',
)
{% endfor %}
{% endif %}

{% endif %}

# COLOUR SETTINGS
ColourSchemeManager(
    name = 'atomtype',
    ranges = ['*', ' C', ' O', ' N', ' S', ' H', ' P'],
    colours = ['grey', 'green', 'red', 'blue', 'yellow', 'grey', 'magenta'],
    colour_wheel_direction = 'clockwise',
)

ParamsManager(
    name = 'gui_params',
    show_axes = 0,
    background_colour = [0, 0, 0],
)
