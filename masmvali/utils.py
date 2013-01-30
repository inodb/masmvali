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
