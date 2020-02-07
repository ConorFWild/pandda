import re

from bamboo.common.command import CommandManager

def phenix_find_tls_groups(pdb_file):
    cmd = CommandManager('phenix.find_tls_groups')
    cmd.add_command_line_arguments(pdb_file)
    #cmd.print_settings()
    ret_code = cmd.run()

    if ret_code != 0:
        print cmd.output
        print cmd.error
        raise Exception('Failed to determine TLS groups: {}'.format(' '.join(cmd.program)))

    regex = re.compile("refinement\.refine\.adp \{([\s\S]*?)\}")
    tls_command = regex.findall(cmd.output)[0]
    tls_selections = [s.strip() for s in tls_command.split('tls =') if s.strip()]

    return tls_selections
