import os, sys
import pkg_resources

try:
    VERSION = pkg_resources.get_distribution("panddas").version
except:
    VERSION = '(developer -- see setup.py file)'

HEADER_TEXT = """
------------------------------------------------------------------>
-         .__               __
-    ____ |__|____    _____/  |_
-   / ___\\|  \\__  \\  /    \\   __\\
-  / /_/  >  |/ __ \\|   |  \\  |
-  \\___  /|__(____  /___|  /__|     A crystallographic toolbox for
-  /_____/         \\/     \\/        protein structure determination
-
------------------------------------------------------------------>

> {program}
{description}
------------------------------------------------------------------>
"""

class module_info:
    name        = 'giant'
    version     = VERSION
    header_text = HEADER_TEXT
