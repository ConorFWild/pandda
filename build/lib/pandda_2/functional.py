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
