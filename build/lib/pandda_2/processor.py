from collections import OrderedDict

import joblib

import luigi
from pandda_2.luigi_sge import SGEJobTask

from libtbx import easy_mp


# class Processor:
#     def __init__(self):
#         pass
#
#     def __call__(self,
#                  func,
#                  args_list,
#                  ):
#
#         results = []
#         for args in args_list:
#             result = func(*args)
#             results.append(result)
#
#         return results
#
#     def repr(self):
#         repr = OrderedDict()
#         return repr

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


class Task(SGEJobTask):
    func = luigi.Parameter()
    output_path = luigi.Parameter()

    def work(self):
        self.func()

    def output(self):
        return luigi.LocalTarget(str(self.output_path))


class ProcessorLuigi:

    def __init__(self,
                 jobs=10,
                 parallel_env="smp",
                 n_cpu=12,
                 run_locally=False,
                 ):
        self.jobs = jobs
        self.parallel_env = parallel_env
        self.n_cpu = n_cpu
        self.run_locally = run_locally

    def __call__(self,
                 funcs,
                 output_paths=None,
                 result_loader=None,
                 shared_tmp_dir=None,
                 ):
        tasks = [Task(func=func,
                      output_path=output_path,
                      shared_tmp_dir="/dls/science/groups/i04-1/conor_dev/pandda/lib-python/pandda/pandda_analyse_dask/luigi_test",
                      parallel_env=self.parallel_env,
                      n_cpu=self.n_cpu,
                      run_locally=False,
                      )
                 for func, output_path
                 in zip(funcs, output_paths)
                 ]

        luigi.build(tasks,
                    local_scheduler=True,
                    workers=self.jobs,
                    detailed_summary=False,
                    )

        if result_loader:
            results = [result_loader(output_path)
                       for output_path
                       in output_paths
                       ]

        else:
            results = []

        return results

    def repr(self):
        repr = OrderedDict()
        repr["jobs"] = self.jobs
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
