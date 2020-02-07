
from bamboo.common.masks import MaskCollection
from bamboo.common import Meta


class HolderList(object):


    _holder_class = None

    def __init__(self):
        """Class for grouping many holders together"""
        self._holder_list = []
        self._masks = None
        self._update_masks()

    def __getitem__(self, idx):
        if   isinstance(idx, int): return self.get(num=idx)
        elif isinstance(idx, str): return self.get(tag=idx)
        else: raise Exception('CANNOT INDEX EXCEPT BY int OR str. TYPE GIVEN: {!s}'.format(type(idx)))

    def __call__(self):
        """Return all holders"""
        return self.all()

    def _get_tag(self, item):
        return item.tag
    def _get_num(self, item):
        return item.num

    def _update_masks(self):
        self._masks = MaskCollection(index_ids=self.all_tags())

    def all(self):
        """Return all holders"""
        return self._holder_list

    def all_nums(self):
        """Return the list of holder ids"""
        return [self._get_num(h) for h in self.all()]

    def all_tags(self):
        """Return the list of holder tags"""
        return [self._get_tag(h) for h in self.all()]

    def all_masks(self):
        """Return the mask object"""
        return self._masks

    def size(self, mask_name=None, invert=False):
        """Return the number of holders in the list (with optional mask applied)"""
        if mask_name: return len(self.mask(mask_name=mask_name, invert=invert))
        else:         return len(self.all())

    def mask(self, mask_name=None, mask=None, invert=False):
        """Retrieve a masked list of datasets"""
        if (mask_name is not None) and (mask is not None):
            raise Exception('Cannot supply both mask_name and mask')
        if mask_name is not None:
            return self._masks.mask_list(name=mask_name, input_list=self.all(), invert=invert)
        elif mask is not None:
            return self._masks.mask_custom(mask=mask,    input_list=self.all(), invert=invert)
        else:
            return self.all()

    def add(self, new):
        """Add new datasets - will reset mask object"""

        for n in new:
            if self._holder_class:
                assert isinstance(n, self._holder_class), 'OBJECTS MUST BE OF TYPE: {!s}\n(ADDED OF TYPE: {!s})'.format(self._holder_class, type(n))
            assert isinstance(self._get_tag(n), str), 'TAG MUST BE str. Type given: {!s}'.format(type(self._get_tag(n)))
            assert self._get_tag(n) not in self.all_tags(), 'HOLDER WITH TAG ALREADY EXISTS: {!s}'.format(self._get_tag(n))
            assert isinstance(self._get_num(n), int), 'NUM MUST BE int. Type given: {!s}'.format(type(self._get_num(n)))
            assert self._get_num(n) not in self.all_nums(), 'HOLDER WITH NUM ALREADY EXISTS: {!s}'.format(self._get_num(n))
            # No problems, add to list
            self._holder_list.append(n)

        self._update_masks()

        return self

    def get(self, tag=None, num=None):
        """Get a dataset by tag or num"""

        assert [num, tag].count(None) == 1, 'Must give EITHER num OR tag'
        if num: matching = [m for m in self.all() if self._get_num(m) == num]
        else:   matching = [m for m in self.all() if self._get_tag(m) == tag]
        if len(matching) == 0: raise Exception('NO MATCHING HOLDER FOUND - NUM: {!s}, TAG: {!s}'.format(num, tag))
        if len(matching) != 1: raise Exception('MORE THAN ONE MATCHING HOLDER FOUND - NUM: {!s}, TAG: {!s}'.format(num, tag))
        return matching[0]

