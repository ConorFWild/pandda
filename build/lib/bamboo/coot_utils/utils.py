
from bamboo.common.command import CommandManager
from bamboo.ccp4_utils import MtzSummary
from bamboo.coot_utils.scripts import *

def validate_coot_script(script):
    """Performs checks on a coot script"""

    contents = open(script, 'r').read().strip()

    assert contents.endswith('coot_real_exit(1)'), 'Coot script does not EXIT! {!s}'.format(script)

    return

def run_coot(script, graphical=False, noguano=True):
    """Runs coot with the provided script"""

    # Check that the script is valid (e.g. causes coot to exit at the end)
    validate_coot_script(script)

    # Stop coot droppings?
    coot_flags = ['--no-guano']*(noguano)
    # Run with/without graphics?
    if graphical:
        coot_flags.extend(['-script',script])
    else:
        coot_flags.extend(['--no-graphics','-s',script])

    # Initialise
    COOT = CommandManager('coot')
    # Load arguments
    COOT.add_command_line_arguments(*coot_flags)
    # Run!
    COOT.run()

    return COOT

