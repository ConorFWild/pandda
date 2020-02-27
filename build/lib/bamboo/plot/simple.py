import numpy

from matplotlib import pyplot
from bamboo.maths import finite

try:
    pyplot.style.use('ggplot')
except:
    pass

def simple_histogram(filename, data, title, x_lab='x', n_bins=30):
    """Quick histogram function"""

    fig = pyplot.figure()
    pyplot.title(title)
    pyplot.hist(x=finite(data), bins=n_bins)
    pyplot.xlabel(x_lab)
    pyplot.ylabel('Count')
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close(fig)

    return

def simple_bar(filename, y_vals, x_labels, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True):
    """Quick bar plot function"""

    assert len(y_vals) == len(x_labels)

    # Plot sequential bars
    x_vals = numpy.arange(len(y_vals))

    fig = pyplot.figure()
    pyplot.title(title)
    pyplot.bar(left=x_vals, height=y_vals, width=1)
    pyplot.xlabel(x_lab)
    pyplot.ylabel(y_lab)
    pyplot.xticks(x_vals+0.5, x_labels)
    pyplot.xlim(x_lim)
    pyplot.ylim(y_lim)
    if rotate_x_labels:
        locs, labels = pyplot.xticks()
        pyplot.setp(labels, rotation=90)
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close(fig)

    return

def simple_boxplot(filename, y_vals, x_labels, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True):

    #assert set(map(len,y_vals)) == {len(x_labels)}

    fig = pyplot.figure()
    pyplot.title(title)
    pyplot.boxplot(y_vals, labels=x_labels, showmeans=True)
    pyplot.xlabel(x_lab)
    pyplot.ylabel(y_lab)
    pyplot.xlim(x_lim)
    pyplot.ylim(y_lim)
    if rotate_x_labels:
        locs, labels = pyplot.xticks()
        pyplot.setp(labels, rotation=90)
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close(fig)

    return

def simple_scatter(filename, x_vals, y_vals, title, x_lab='x', y_lab='y', x_lim=None, y_lim=None, rotate_x_labels=True):
    """Quick scatter plot function"""

    assert len(x_vals) == len(y_vals)

    fig = pyplot.figure()
    pyplot.title(title)
    pyplot.scatter(x=x_vals, y=y_vals)
    pyplot.xlabel(x_lab)
    pyplot.ylabel(y_lab)
    pyplot.xlim(x_lim)
    pyplot.ylim(y_lim)
    if rotate_x_labels:
        locs, labels = pyplot.xticks()
        pyplot.setp(labels, rotation=90)
    pyplot.tight_layout()
    pyplot.savefig(filename)
    pyplot.close(fig)

    return

