class WrapDebug:
    def __init__(self,
                 func,
                 debug=True,
                 ):
        self.func = func
        self.debug = debug

    def __call__(self, *args, **kwargs,):

        if self.debug:
            print("Starting function: {}".format(self.func))

        results = self.func(*args, **kwargs)

        if self.debug:
            print("Finished function: {}".format())

        return results