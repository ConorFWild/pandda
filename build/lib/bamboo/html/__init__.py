
import jinja2
import urlparse, urllib

BAMBOO_HTML_ENV = jinja2.Environment(loader=jinja2.PackageLoader('bamboo', 'resources'))

def path2url(path):
    return urlparse.urljoin('file:', urllib.pathname2url(path))

def png2base64str(path):
    """Convert a png for embedding into an html page"""
    contents = open(path, 'rb').read().encode('base64').replace('\n', '')
    return contents
