
settings_phil = """
settings
    .help = "General Settings"
{
    cpus = 1
        .type = int
    verbose = True
        .type = bool
    plot_graphs = True
        .help = "Output graphs using matplotlib"
        .type = bool
    plotting {
        backend = 'agg'
            .help = "Backend to use in matplotlib"
            .type = str
    }
}
"""

image_phil = """
image
    .help = "Settings for image generation"
{
    size = 300x200
        .help = "image size of width x height (inches)"
        .type = str
    format = *png jpg svg
        .help = "select output image format
        .type = choice
}
"""
