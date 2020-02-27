import operator
import numpy, pandas


class MaskCollection(object):


    def __init__(self, index_ids=None):
        """Stores lots of similar mask objects. Masks can be combined in multiple ways"""
        self.table = pandas.DataFrame()
        self._index_length = None
        if index_ids is not None:
            self.set_index_ids(index_ids)

    def invert_mask(self, mask):
        """Invert a list of boolean values"""
        mask = numpy.array(mask).astype(bool)
        assert mask.ndim == 1
        mask = map(operator.not_, mask)
        return pandas.Series(data=mask, index=self.table.index, dtype=bool)

    def set_index_length(self, length):
        """Set the mask length"""
        assert not self._index_length
        self._index_length = length

    def set_index_ids(self, ids):
        assert not self.table.index.name
        if not self._index_length:
            self.set_index_length(len(ids))
        assert len(ids) == self._index_length
        self.table['id'] = ids
        self.table.set_index('id', inplace=True)

    def add_mask(self, name, values, overwrite=False):
        """Add a new mask. Values can be list of bools or bool"""
        if not isinstance(values, bool):
            # Convert to numpy array
            values = numpy.array(values).astype(bool)
            # Allow to be a "2D" column array
            if values.shape == (self._index_length, 1):
                values = values.flatten()
            assert values.ndim == 1
            assert len(values) == self._index_length, 'if values is list must have same length as index_length'
        if not overwrite:
            assert name not in self.table.columns, 'to overwrite columns require overwrite=True'
        self.table[name] = values

    def get_mask(self, name, invert=False):
        """Get a single mask"""
        mask = self.table[name]
        if invert: return self.invert_mask(mask)
        else:      return mask

    def has_mask(self, name):
        """Check if it has a single mask"""
        return name in self.table.columns

    def set_value(self, name, id, value):
        """Set a particular entry in the mask `mask_name`, corresponding to `id` to `value`"""
        assert name in self.table.columns
        assert id in self.table.index
        self.table[name][id] = bool(value)

    def get_value(self, name, id):
        """Set a particular entry in the mask `mask_name`, corresponding to `entry_id` to `value`"""
        assert name in self.table.columns
        assert id in self.table.index
        return self.table[name][id]

    def mask_index(self, name, invert=False):
        """Return the index ids selected by a mask"""
        return self.mask_custom(mask=self.table[name], input_list=self.table.index, invert=invert)

    def mask_list(self, name, input_list, invert=False):
        """Mask input_list according to the mask with name mask_name. Invert mask with `invert`"""
        return self.mask_custom(mask=self.table[name], input_list=input_list, invert=invert)

    def mask_custom(self, mask, input_list, invert=False):
        """Mask input_list according to the mask. Invert mask with `invert`"""
        assert isinstance(invert, bool)
        assert len(mask) == self._index_length, 'Given mask is the wrong length'
        assert len(input_list) == self._index_length, 'Given list is the wrong length'
        return [l for m, l in zip(mask, input_list) if m != invert]

    def combine_masks(self, names, invert=False, operation='or', invert_output=False):
        """Combine multiple masks. Invert each mask by invert (can be list) and then join using booling operation 'and' or 'or'"""
        assert operation in ['or','and']
        # Expand invert to list if given as bool
        if isinstance(invert, bool):
            invert = [invert]*len(names)
        else:
            assert len(invert) == len(names)
        # Combines masks into one table
        masks = pandas.concat([self.get_mask(n, invert=i) for n,i in zip(names,invert)], axis=1)
        # Select the operation: "and" is "min" since all must be true for it to be true
        if operation == 'and':
            func = min
        else:
            func = max
        # Combine the masks using and or or
        out_mask = masks.apply(func, axis=1)
        if invert_output: out_mask = self.invert_mask(out_mask)
        return out_mask

