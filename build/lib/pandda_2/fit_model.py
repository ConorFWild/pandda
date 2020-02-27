from collections import OrderedDict

class FitModel:
    def __init__(self):
        pass

    def __call__(self, model, train, test):
        fit_model = model.fit(train, test)

        return fit_model

    def repr(self):
        repr = OrderedDict()
        return repr