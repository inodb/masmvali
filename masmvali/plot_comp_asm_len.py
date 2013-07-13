#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse

import assembly
from plot import savefig_multiple_ext


def plot_length_roc_style(assemblies, names, colors, linestyles, linewidths, cut_off=100):
    ax = plt.subplot(1, 1, 1)
    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        nr_of_contigs = len(cla[cla > cut_off])
        # x contig lengths, 100 dots
        xarray = np.arange(nr_of_contigs, step=nr_of_contigs / 100)
        # sum of bases
        yarray = np.cumsum(cla[cla > cut_off])[xarray]
        plt.plot(xarray, yarray, color=colors[i], label=names[i], linestyle=linestyles[i], linewidth=linewidths[i])
    # Labels
    plt.xlabel('Number of contigs from biggest to smallest')
    plt.ylabel('Sum of bases')
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='lower right')


def plot_length_cum(assemblies, names, colors, shapes, linewidths, cut_off=100):
    for i, asm in enumerate(assemblies):
        contigs = assembly.ContigDict()
        contigs.parse_fasta_lengths(asm)
        contiglens = [c.length for c in contigs.itervalues()]
        contiglens.sort()
        contiglens.reverse()
        cla = np.array(contiglens)
        plt.plot(cla[cla > cut_off], np.cumsum(cla[cla > cut_off]), shapes[i], color=colors[i], label=names[i], linewidth=linewidths[i])
    plt.gca().invert_xaxis()


def main(outputbase, assemblies, names, colors, shapes, linewidths, cut_off=100, roc=False):
    # Change figure size
    plt.rcParams['figure.figsize'] = 8, 6
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 100

    if not roc:
        plot_length_cum(assemblies, names, colors, shapes, linewidths, cut_off)
    else:
        plot_length_roc_style(assemblies, names, colors, shapes, linewidths, cut_off)

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
                               help='Label names of the assembly. Input should be a file with one line for each assembly')
    parser.add_argument('--colors', default=False,
                               help='Plot with given colors, use matplotlib color names e.g. "r" for red. Input should be a file with one line for each assembly')
    parser.add_argument('--shapes', default=False,
                               help='Plot names with given shapes, use matplotlib shape names e.g. "." for point. Input should be a file with one line for each assembly')
    parser.add_argument('--linewidths', default=False,
                               help='Plot names with given linewidth. Input should be a file with one line for each assembly')
    parser.add_argument('--roc', action='store_true',
                               help='Plot roc style')

    args = parser.parse_args()
    if not args.names or not args.colors or not args.shapes:
        raise(Exception('Names, colors and shapes option are required options'))

    args.names = [l.rstrip('\n') for l in open(args.names).readlines()]
    if not len(args.assemblies) == len(args.names):
        raise(Exception('Number of names should be equal to amount of fasta files given.'))

    args.colors = [l.rstrip('\n') for l in open(args.colors).readlines()]
    if not len(args.assemblies) == len(args.colors):
        raise(Exception('Number of colors should be equal to amount of fasta files given.'))

    args.shapes = [l.rstrip('\n') for l in open(args.shapes).readlines()]
    if not len(args.assemblies) == len(args.shapes):
        raise(Exception('Number of shapes should be equal to amount of fasta files given.'))

    if not args.linewidths:
        args.linewidths = [1.0] * len(args.names)
    else:
        args.linewidths = [float(l.rstrip('\n')) for l in open(args.linewidths).readlines()]
        if not len(args.assemblies) == len(args.linewidths):
            raise(Exception('Number of linewidths should be equal to amount of fasta files given.'))

    main(args.outputbase, args.assemblies, args.names, args.colors, args.shapes, args.linewidths, roc=args.roc)
