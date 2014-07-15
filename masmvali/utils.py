import os
import errno

def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def print_dict2tsv(d, filepath):
    with open(filepath, "w") as fh:
        fh.write("\t".join(d.keys()) + "\n")
        fh.write("\t".join(str(v) for v in d.itervalues()) + "\n")


def readtable(tablefile, sep=None):
    """Reads given table separated 'sep'.

    Returns a two-dimensional dictionary. Outer dictionary uses first column in
    'tablefile' as key, inner dictionary uses column names found in the first
    line as column names. The number of column names should be equal to number
    of columns - 1."""
    table2d = {}

    # Get column names
    tfh = open(tablefile, "r")
    line = tfh.readline()
    cols = line.strip().split(sep)

    # Insert rows
    for line in tfh:
        splits = line.strip().split(sep)
        table2d[splits[0]] = {}
        for i in range(1, len(splits)):
            table2d[splits[0]][cols[i - 1]] = splits[i]

    return(table2d)


def read_2col_table(tablefile, sep=None):
    """Reads given table separated 'sep'.

    Return a dictionary where first column become keys and the second
    values."""
    table = {}

    # Get column names
    tfh = open(tablefile, "r")

    # Insert rows
    for line in tfh:
        splits = line.strip().split(sep)
        if len(splits) != 2:
            raise(Exception("Expected table with two columns. Line looks like:\n" + line))
        table[splits[0]] = splits[1]

    return(table)
