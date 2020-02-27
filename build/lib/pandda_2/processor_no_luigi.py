from collections import OrderedDict

import joblib

from libtbx import easy_mp



class Processor:
    def __init__(self):
        pass

    def __call__(self,
                 funcs,
                 output_paths=None,
                 result_loader=None,
                 shared_tmp_dir=None,
                 ):
        results = []
        for func in funcs:
            result = func()
            results.append(result)

        return results

    def repr(self):
        repr = OrderedDict()
        return repr


class ProcessorDict:
    def __init__(self):
        pass

    def __call__(self,
                 funcs,
                 ):
        results = {}
        for key, func in funcs.items():
            result = func()
            results[key] = result

        return results

    def repr(self):
        repr = OrderedDict()
        return repr


class ProcessorDictJoblib:
    def __init__(self,
                 cpus=21,
                 verbosity=8,
                 ):
        self.cpus = cpus
        self.verbosity = verbosity

    def __call__(self,
                 funcs,
                 ):
        keys = funcs.keys()
        values = [funcs[key] for key in keys]

        results = joblib.Parallel(n_jobs=self.cpus,
                                  verbose=self.verbosity,
                                  )(joblib.delayed(value)()
                                    for value
                                    in values
                                    )

        results_dict = {key: results[i]
                        for i, key
                        in enumerate(keys)
                        }

        return results_dict

    def repr(self):
        repr = OrderedDict()
        repr["cpus"] = self.cpus
        repr["verbosity"] = self.verbosity
        return repr


class ProcessorJoblib:
    def __init__(self,
                 cpus=21,
                 verbosity=8,
                 ):
        self.cpus = cpus
        self.verbosity = verbosity

    def __call__(self, funcs):
        results = joblib.Parallel(n_jobs=self.cpus,
                                  verbose=self.verbosity,
                                  )(joblib.delayed(func)()
                                    for func
                                    in funcs
                                    )

        return results


def wrap_call(f):
    return f()


class ProcessorDictEasyMP:
    def __init__(self,
                 cpus=21,
                 verbosity=8,
                 ):
        self.cpus = cpus
        self.verbosity = verbosity

    def __call__(self,
                 funcs,
                 ):
        keys = funcs.keys()
        values = [funcs[key] for key in keys]

        results = easy_mp.pool_map(fixed_func=wrap_call,
                                   args=values,
                                   processes=int(self.cpus),
                                   )

        results_dict = {key: results[i]
                        for i, key
                        in enumerate(keys)
                        }

        return results_dict

    def repr(self):
        repr = OrderedDict()
        repr["cpus"] = self.cpus
        repr["verbosity"] = self.verbosity
        return repr
