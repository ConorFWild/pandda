from collections import OrderedDict

class EvaluateModel:
    def __init__(self):
        pass

    def __call__(self, model, xmap,):

        zmap = model.evaluate(xmap)

        return zmap

    def repr(self):
        repr = OrderedDict()
        return repr