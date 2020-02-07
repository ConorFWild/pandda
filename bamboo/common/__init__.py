
class Report(object):
    """Used to compile message strigns from functions -- useful for returning function history upon errors"""

    def __init__(self, s='', verbose=False):
        self.s = s
        self.bar = '-'*100
        self.verbose = verbose
        if self.s and self.verbose: print self.s

    def __call__(self, s):
        if self.verbose: print s
        self.s += '\n' + s

    def __str__(self):
        return self.bar + '\n' + self.s + '\n' + self.bar


class ListStream:
    """Replaces a file object for functions that print to screen - writes to list of lines instead"""


    def __init__(self):
        self.data = []

    def __str__(self):
        return ''.join(self.data)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        return iter(self.data)

    def write(self, s):
        self.data.append(s)


class Meta(object):
    """Object for storing random data - can be edited and added to as required"""


    _initialized = False

    def __init__(self, args=None):
        if isinstance(args, dict):
            for k in args:
                self.__dict__[k] = args[k]
        elif isinstance(args, list):
            for l in args:
                self.__dict__[l] = None
        elif args is not None:
            raise Exception('args must be dict or list')
        _initialized = True

    def summary(self):
        out = []
        for k in self.__dict__:
            out.append('{!s}: {!s}'.format(k, self.__dict__[k]))
        return '\n'.join(out)


class Info(Meta):
    """Same as Meta, but cannot change variables after setting _initialized=True"""


    def __setattr__(self, name, value):
        if self._initialized and (not hasattr(self, name)):
            raise AttributeError('Cannot set new attributes after initialisation: {!s}'.format(name))
        object.__setattr__(self, name, value)



