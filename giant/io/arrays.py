import os, sys
import numpy

from bamboo.common.status import status_bar

def write_numpy_array_to_file(np_array, output_file, status=True):
    """Takes numpy array and writes to file. Chooses axis of numpy array so as to maximise the number of rows of the output file"""

    assert isinstance(np_array, numpy.array), 'Array must be numpy array!'
    assert not os.path.exists(output_file), 'Output file already exists!'
    assert np_array.ndim <= 2, 'Dimension of Matrix must be less than 3!'

    # Find the longest dimension of the array
    long_num = max(np_array.shape)
    long_dim = np_array.shape.index(row_num)

    # Array is wider than it is long - need to transpose
    if long_dim == 1:
        np_array=np_array.swapaxes(0,1)

    # Now iterate through and write out the lines
    with open(output_file, 'w') as fh:
        for row in xrange(long_num):
            if status: status_bar(n=row, n_max=long_num)
            output_line = [str(round(x,3)) for x in np_array[row]]
            fh.write(', '.join(output_line)+'\n')
