from argparse import ArgumentParser

class ArgParser(ArgumentParser):
    """Minor Modifications to Argument Parser"""
    def parse_arguments(self, args=None):
        """Parse Arguments and return results also as a dictionary"""
        if args:
            self.inputs = self.parse_args(args)
        else:
            self.inputs = self.parse_args()

        return self.inputs, self.inputs.__dict__

