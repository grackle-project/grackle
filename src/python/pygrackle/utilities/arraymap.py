########################################################################
#
# Implements the arraymap class
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import itertools
import numpy as np

def _memory_region_range(arr):
    # returns a 2-tuple specifying the range of the memory addresses (i.e. the
    # values of a c-pointer) that bound the memory region used by `arr`
    #
    # The first element specifies the minimum address of any element in arr
    # The holds the minimum memory address outside of the memory region (it is
    # larger than the first element)
    ptr_first_elem = arr.__array_interface__['data'][0]
    strides, shape, itemsize = arr.strides, arr.shape, arr.itemsize

    if len(strides) == 1 and strides[0] >= 0:
        ptr_last_elem = ptr_first_elem + strides[0] * (shape[0] - 1)
        return (ptr_first_elem, ptr_last_elem + itemsize)
    
    # the general case is messy (strides can have mix of positives & negatives)
    index_offsets = [] # offsets from the address of very first index
    for index in itertools.product(*[(0, axlen - 1) for axlen in shape]):
        index_offsets.append(int(np.dot(index, strides)))
    return (ptr_first_elem + min(index_offsets),
            ptr_first_elem + max(index_offsets) + arr.itemsize)


class ArrayMap:
    """
    An arraymap object represents a dict of arrays, which share the same dtype,
    shape, and layout in memory.

    Unlike a normal dict, keys can't be added/removed after initialization.

    Parameters
    ----------
    obj : mapping
        A mapping of numpy arrays. Each contained array must share a common
        shape. We could make this more flexible and allow the mapping to hold
        arraylike objects instead.
    dtype : data-type, optional
        The shared datatype of the contained arrays
    copy : bool, optional
        Specifies whether to make deepcopies of all entries in obj. In a
        departure from numpy behavior, an error is raised when this is `False`
        and an arraymap cannot be initialized with shallowcopies. When this is
        `None`, (the default value), efforts are made to use shallowcopies, but
        deepcopies are used if necessary.
    """
    def __init__(self, obj, *, dtype = "f8", copy = None):

        # stage 1: prepare all variables that will be used or modified in
        #          subsequent parts of this function
        dtype = np.dtype(dtype)
        shape = None
        # the following will be assigned an explanatory error message if there
        # is some reason we can't reuse existing arrays (i.e. if there's some
        # reason we NEED to make deepcopies)
        _reuse_err = None

        # stage 2: Perform some checks on obj and query some information
        if isinstance(obj, ArrayMap):
            # bypass most checks -- obj must already satisfy most requirements
            if dtype != obj.dtype:
                _reuse_err = "obj.dtype isn't the same as the specified dtype"
            shape = obj.array_shape
        else: # general case
            # for now, we are intentionally overly restrictive (it's easy to
            # become more permissive later)

            # _address_ranges is filled with 3-tuples (ptr_min, ptr_stop, key)
            #  - ptr_min holds the lowest memory address of an array element
            #  - ptr_stop holds the first memory address after ptr_min not used
            #    by an array
            #  - key records the array in question
            _address_ranges = []
            if (copy is not None) and not bool(copy):
                def perform_reuse_check(key, arr):
                    start_addr, stop_addr = _memory_region_range(arr)
                    _address_ranges.append((start_addr, stop_addr))
                    if not arr.flags['C_CONTIGUOUS']:
                        return f"array associated with {key!r} not contiguous"
                    elif dtype != arr.dtype:
                        return (f"array associated with {key!r} doesn't have "
                                "specified dtype")
            else:
                perform_reuse_check = lambda key, arr : None

            for key, arr in obj.items():  # make an error-checking pass
                if shape is None:
                    shape = arr.shape

                if type(arr) is not np.ndarray:
                    # we  use isinstance to make sure we aren't
                    # using a subclass like unyt.unyt_array (maybe allow that
                    # in future)
                    raise ValueError("all values of obj must be instances of "
                                     "np.ndarray (subclasses are not allowed)")
                elif shape != arr.shape:
                    raise ValueError(f"arr associated with {key!r} has an "
                                     "unexpected shape")
                elif _reuse_err is None:
                    _reuse_err = perform_reuse_check(key, arr)
                    # do NOT break out of loop early (for `copy is None` case)

            # perform extra check for whether we can directly reuse arrays
            if (_reuse_err is None) and (len(_address_ranges) > 1):
                _address_ranges.sort()
                last_stop = 0
                for start, stop in _address_ranges:
                    if start < last_stop:
                        _reuse_err = "there is overlap between arrays"
                    last_stop = stop

        # stage 3: Use queried info to perform some final checks and, if
        #          applicable, make final judgement on directly reusing arrays
        if shape is None:
            raise ValueError("When obj is empty, the shape arg is required")
        elif len(shape) != 1:
            raise NotImplementedError("There's currently no support for "
                                      "arrays that aren't 1D")
        if copy is None:
            copy = _reuse_err is None
        elif (not bool(copy)) and _reuse_err is not None:
            raise ValueError("Can't reuse the arrays from obj argument to "
                             "construct arraymap: " + _reuse_err)

        # stage 4: finally, initialize self
        if copy:
            self._data = dict((k, np.array(a, dtype, copy = True, order = 'C'))
                              for k,a in obj.items())
        else:
            self._data = dict((k,a) for (k,a) in obj.items())
        self._dtype = dtype
        self._array_shape = shape

    @property
    def dtype(self):
        """Data-type of elements of contained arrays
        """
        return self._dtype

    @property
    def array_shape(self):
        """Shape of each contained array
        """
        return self._array_shape

    @classmethod
    def fromkeys(cls, keys, *, shape, dtype = "f8", fill_value = None):
        """
        Construct a new arraymap with keys from the `keys` argument and values
        set to arrays with properties specified by the other kwargs

        Notes
        -----
        The name is inspired by `dict.fromkeys`
        """
        # todo: consider coming up with a more efficient implementation
        #       (i.e. how to bypass unnecessary tests in the constructor)

        if fill_value is None:
            fn = lambda : np.empty(shape, dtype = dtype)
        else:
            fn = lambda : np.full(shape, fill_value, dtype = dtype)

        return cls( dict((k,fn()) for k in keys) )

    def __getitem__(self, key):
        """Returns the array associated with `key`
        """
        return self._data.__getitem__(key)

    def __contains__(self, key):
        """Returns whether `self` contains the specified `key`
        """
        return self._data.__contains__(key)

    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return self._data.__len__()

    def keys(self):
        """Returns a view of the keys
        """
        return self._data.keys()

    def items(self):
        return self._data.items()

    def values(self):
        return self._data.values()

    def get(self, *args):
        return self._data.get(*args)

    def __eq__(self, other):
        raise NotImplementedError(
            "We have chosen not to implement this method. It's unclear whether "
            "we should compare contents or identities of self.values() and "
            "other.values()")

    def __ne__(self, other):
        return not self.__eq__(other)
