from functools import reduce


def chain(data, funcs):

    funcs.insert(0, data)

    result = reduce(lambda d, f: f(d), funcs)

    return result


def curry_kwargs(func, **kwargs):
    return lambda *args: func(*args, **kwargs)


def multimap(func, *args):

    func_to_map = lambda x: func(*x)

    results = map(func_to_map,
                  zip(args),
                  )

    return results


def mapdict(func, dict):
    items = dict.items()
    keys = [item[0] for item in items]
    values = [item[1] for item in items]
    results = map(func, values)
    results_dict = dict(zip(keys, results))
    return results_dict


def mapdictparallel(func, dictionary):
    keys = dictionary.keys()
    values = dictionary.values()

    results = joblib.Parallel(n_jobs=21, verbose=10, batch_size=1)(
        joblib.delayed(func)(value) for value in values
    )

    results_dict = dict(zip(keys, results))
    return results_dict
