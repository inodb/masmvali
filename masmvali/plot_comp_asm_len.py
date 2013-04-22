#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse

import assembly
from plot import savefig_multiple_ext


def main(outputbase, assemblies, names, colors, shapes, cut_off=100):
    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        plt.plot(cla[cla > cut_off], np.cumsum(cla[cla > cut_off]), shapes[i], color=colors[i], label=names[i])
    plt.gca().invert_xaxis()

    savefig_multiple_ext(outputbase)

    # Return png name of outputted file
    return outputbase + ".png"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('outputbase', metavar='outputbase',
                               help='Output file basename, outputted in png, svg and jpg')
    parser.add_argument('assemblies', metavar='assembly.fa',
                        nargs='+',  help='Contigs or scaffolds from an assembly')
    parser.add_argument('--names', default=False,
                               help='Label given assembly files by given names, separate with ,')
    parser.add_argument('--colors', default=False,
                               help='Plot names with given colors, use matplotlib color names e.g. "r" for red, separate multiple by ",", should be equally long as names')
    parser.add_argument('--shapes', default=False,
                               help='Plot names with given shapes, use matplotlib shape names e.g. "." for point, separate multiple by ",", should be equally long as names')

    args = parser.parse_args()
    if not args.names or not args.colors or not args.shapes:
        raise(Exception('Names, colors and shapes option are required options'))
    args.names = args.names.split(",")
    if not len(args.assemblies) == len(args.names):
        raise(Exception('Number of names should be equal to amount of fasta files given.'))
    args.colors = args.colors.split(",")
    if not len(args.assemblies) == len(args.colors):
        raise(Exception('Number of colors should be equal to amount of fasta files given.'))
    args.shapes = args.shapes.split(",")
    if not len(args.assemblies) == len(args.shapes):
        raise(Exception('Number of shapes should be equal to amount of fasta files given.'))
    main(args.outputbase, args.assemblies, args.names, args.colors, args.shapes)
