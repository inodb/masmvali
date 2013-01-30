#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

from plot import savefig_multiple_ext
from utils import make_dir

MISSING_VALUE = "N/A"
DELIMITER = "\t"


def plot_per_row_legend_n_save(a, x, y, namefield, xlab, ylab, outdir, baseout=False, xlim=False, ylim=False, size_field=False):
    """Plot two columns of a numpy array against each other. Give each point a
    different color and include a legend."""
    # Plot
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    names = np.unique(a[namefield]).tolist()
    assert len(names) < len(colors)
    colord = dict(zip(names, colors[:len(names)]))
    ax = plt.subplot(1, 1, 1)
    if not size_field:
        for i in range(len(a)):
            ax.plot(a[x][i], a[y][i], '.', color=colord[a[namefield][i]], label=a[namefield][i], markersize=20)
    else:
        for i in range(len(a)):
            ax.plot(a[x][i], a[y][i], '.', color=colord[a[namefield][i]], label=a[namefield][i], markersize=int(a[size_field][i] / float(max(a[size_field])) * 20))
    handles, labels = ax.get_legend_handles_labels()
    labelsu = list(set(labels))
    handlesu = [handles[labels.index(l)] for l in labelsu]
    ax.legend(handlesu, labelsu, numpoints=1, loc=0)

    if xlim:
        plt.xlim(xmin=xlim[0], xmax=xlim[1])
    else:
        # X and Y are plotted on the direct ends increase plot size
        xmin, xmax = plt.xlim()
        xs = (xmax - xmin) / 10
        assert xs > 0
        plt.xlim(xmin - xs, xmax + xs)

    if ylim:
        plt.ylim(ymin=ylim[0], ymax=ylim[1])

    plt.xlabel(xlab)
    plt.ylabel(ylab)

    # Save
    if not baseout:
        baseout = "{x}_vs_{y}".format(x=x, y=y)

    savefig_multiple_ext(outdir + "/" + baseout)

    return baseout + ".png"


def main(asmstats, outdir, size_field=False, names=False):
    #DTYPE = [('sum_purest_bases', '<f8'), ('sum_bases', '<i8'), ('kmer_type', '|S3'), ('n50', '<i8'), ('trim_n', '<i8'), ('kmax', '|S3'), ('trim_n_mapping', '<i8'), ('l50', '<i8'), ('name', '|S6'), ('global_purity', '<f8'), ('cut_off', '<i8'), ('aln_purity', '<f8'), ('kmin', '|S3'), ('sum_ref_lengths', '<i8'), ('metagenome_cov', '<f8'), ('kmer_size', '|S3'), ('aln_ratio', '<f8'), ('asm_type', '|S3'), ('max_contig_length', '<i8')]
    a = np.genfromtxt(asmstats, names=True, dtype=None, delimiter=DELIMITER, missing_values=MISSING_VALUE, usemask=True)

    #avalues = []
    #for asms in asmstats:
    #    avalues.append(tuple(open(asms).readlines()[1].strip().split(DELIMITER)))
    #a = np.array(avalues, dtype=DTYPE)

    # If names are specified, filter by given names
    if names:
        a = a[np.in1d(a["name"], names)]

    make_dir(outdir)
    p1 = plot_per_row_legend_n_save(a, "l50", "global_purity", "name", "L50", "Global purity", outdir, ylim=[0, 1], size_field=size_field)
    p2 = plot_per_row_legend_n_save(a, "l50", "aln_purity", "name", "L50", "Alignment purity", outdir, ylim=[0, 1], size_field=size_field)
    p3 = plot_per_row_legend_n_save(a, "l50", "aln_ratio", "name", "L50", "Alignment ratio", outdir, ylim=[0, 1], size_field=size_field)
    p4 = plot_per_row_legend_n_save(a, "l50", "metagenome_cov", "name", "L50", "Metagenome coverage", outdir, ylim=[0, 1], size_field=size_field)

    # Output plots in HTML
    sdir = os.path.dirname(os.path.realpath(__file__))
    template = open(sdir + '/template/cmp-asm-template.html').read()
    with open(outdir + '/index.html', 'w') as fh:
        fh.write(template.format(assemblies=" ".join(np.unique(a["name"]).tolist()),
                                 plot1=p1,
                                 plot2=p2,
                                 plot3=p3,
                                 plot4=p4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('asmstats', metavar='asm-stats.tsv',
                               help='asm-stats.tsv output from validation.py')
    parser.add_argument('outdir', metavar='outdir',
                               help='Output directory')
    parser.add_argument('--sizefield', default=False,
                               help='Scale points according to size in given column')
    parser.add_argument('--names', default=False,
                               help='Only plot rows with given names, seperate multiple by ","')
    args = parser.parse_args()
    out = args.outdir.rstrip('/')
    if args.names:
        args.names = args.names.split(",")
        if not len(args.names) > 1:
            raise(Exception('Less than 1 name given for --names? Why compare a single assembly?'))
    main(args.asmstats, out, args.sizefield, args.names)
