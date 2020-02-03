
import jinja2

# TODO Change and fix

# PANDDA_HTML_ENV = jinja2.Environment(loader=jinja2.PackageLoader('pandda', 'templates'))

# PANDDA_HTML_ENV = jinja2.Environment(loader=jinja2.FileSystemLoader('/dls/science/groups/i04-1/conor_dev/pandda/lib-python/pandda/templates'))

# PANDDA_HTML_ENV = jinja2.Environment(loader=jinja2.FileSystemLoader('/dls/science/groups/i04-1/conor_dev/pandda/lib-python/pandda/templates'))

PANDDA_HTML_ENV = jinja2.Environment(loader=jinja2.PackageLoader('pandda_inspect', 'templates'))

