from subprocess import Popen, PIPE
import sys
import os
import errno


def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def get_shell_output(command, silent=True):
    """Run command using Popen. Return a tuple of the output, error and the
       return code of given command."""
    if not silent:
        print(command)
    event = Popen(command, shell=True, stdin=PIPE, stdout=PIPE)
    output, error = event.communicate()
    if not silent:
        sys.stdout.write(output)
    return output, error, event.returncode

def get_outdir():
    return os.path.dirname(os.path.realpath(__file__)) + '/'

def get_testfile(name):
    sdir = os.path.dirname(os.path.realpath(__file__)) + '/'

    return os.path.abspath(sdir + {
        'cm-500pgun-asm-b2mv31-fa' : 'data/chris-mock/Sample_500pg_unbalanced/ma-out/assemblies/metavelvet/noscaf_31/bambus2/bambus2.scaffold.linear.fasta',
        'cm-500pgun-asm-b2mv31-bam': 'data/chris-mock/Sample_500pg_unbalanced/ma-out/assemblies/metavelvet/noscaf_31/bambus2/val/map/bambus2.scaffold.linear.fasta_reads-smds.bam',
        'cm-500pgun-val-nucmer'    : 'data/chris-mock/Sample_500pg_unbalanced/ma-out/assemblies/metavelvet/noscaf_31/bambus2/val/nucmer.coords',
        'cm-500pgun-ref-stats'     : 'data/chris-mock/Sample_500pg_unbalanced/ma-out/reference-stats/ref.stats',
        'cm-500pgun-ref-bam'       : 'data/chris-mock/Sample_500pg_unbalanced/ma-out/reference-stats/ref_500pg_unbalanced.qtrim-smds.bam',
        'cm-ref-phyl'              : 'data/chris-mock/reference/phylogeny-references.tsv',
        'cm-ref-fa'                : 'data/chris-mock/reference/references.fa'
    }[name])
