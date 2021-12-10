import param


class Cache(param.Parameterized):
    pass


class MemoryCache(param.Parameterized):

    _cache = param.Dict(default={})

    def __getitem__(self, item):
        return self._cache.get(item)

    def __setitem__(self, key, value):
        self._cache[key] = value
