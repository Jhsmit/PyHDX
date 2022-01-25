import param
import pandas as pd


class Cache(param.Parameterized):
    def __getitem__(self, item):
        return None

    def __setitem__(self, key, value):
        pass

    def __contains__(self, item):
        return False


class MemoryCache(Cache):

    _cache = param.Dict(default={})

    max_items = param.Integer(None, doc="Maximum number of items allowed in the cache")

    def __getitem__(self, item):
        return self._cache.__getitem__(item)

    def __setitem__(self, key, value):
        if self.max_items is not None and len(self._cache) >= self.max_items:
            self._cache.popitem()

        self._cache[key] = value

    def __contains__(self, item):
        return item in self._cache


class HybridHDFCache(Cache):
    """

    Hybrid HDFStore / Memory cache

    Sometimes there are errors depending on the dtypes of dataframes stored

    """

    file_path = param.String()

    _store = param.ClassSelector(class_=pd.HDFStore)

    _cache = param.Dict(default={})

    bytes_threshold = param.Integer(default=int(1e8))

    def __init__(self, **params):
        super().__init__(**params)
        if self.file_path is not None:
            self._store = pd.HDFStore(self.file_path)

    def __getitem__(self, item):
        key = str(item)
        try:
            return self._cache.__getitem__(key)
        except KeyError:
            return self._store.__getitem__(key)

    def _store_put(self, key, value):
        try:
            self._store[key] = value

            # Check if reading back the dataframe works
            try:
                _value = self._store[key]
            except AttributeError:
                del self._store[key]
                self._cache[key] = value

        except (
            NotImplementedError,
            TypeError,
        ):  # pytables does not support categorical dtypes
            self._cache[key] = value

    def __setitem__(self, key, value):
        key = str(key)
        if (
            isinstance(value, pd.DataFrame)
            and value.memory_usage().sum() > self.bytes_threshold
        ):
            self._store_put(key, value)
        elif (
            isinstance(value, pd.Series) and value.memory_usage() > self.bytes_threshold
        ):
            self._store_put(key, value)
        else:
            self._cache[str(key)] = value

    def __contains__(self, item):
        return str(item) in self._cache.keys() | self._store.keys()

    # todo with statement for creating caches?
    # def __exit__(self):
    #     pass
