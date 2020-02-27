
class bcolours:
    ENDC = '\033[0m'
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class pretty_string(object):
    def __init__(self, string):
        """Extended class for formatting strings"""
        self.string = string
    def __repr__(self):
        return str(self.__class__)+':'+ self.__str__()
    def __str__(self):
        return self.string
    def __add__(self, other):
        return pretty_string(self.string+other.string)

    def _apply_flag(self, flag):
        """Add flag to the string, inserting after all other ENDC strings"""
        string = self.string[:self.string.rfind(bcolours.ENDC)] if self.string.endswith(bcolours.ENDC) else self.string
        return pretty_string(flag+string.replace(bcolours.ENDC,bcolours.ENDC+flag)+bcolours.ENDC)

    def header(self):
        return self._apply_flag(bcolours.HEADER)
    def bold(self):
        return self._apply_flag(bcolours.BOLD)
    def underline(self):
        return self._apply_flag(bcolours.UNDERLINE)

    def warning(self):
        return self._apply_flag(bcolours.WARNING)
    def fail(self):
        return self._apply_flag(bcolours.FAIL)

    def red(self):
        return self.fail()
    def yellow(self):
        return self.warning()
    def green(self):
        return self._apply_flag(bcolours.OKGREEN)
    def blue(self):
        return self._apply_flag(bcolours.OKBLUE)
    def purpink(self):
        return self._apply_flag(bcolours.HEADER)




