#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
from utils import make_dir


MISSING_VALUE = "N/A"
DELIMITER = "\t"


def savefig_multiple_ext(baseout):
    """Save figure with multiple extensions and close afterwards."""
    plt.savefig(baseout + ".pdf")
    plt.savefig(baseout + ".svg")
    plt.savefig(baseout + ".png")
    plt.close()


def percent_barplot(barvalues, outdir='.', title="Barplot", baseout="barplot", vlab='#'):
    """Plot a list of tuples [(barname, value), ...]" in a barplot. All values
    are expressed as a percentage of the sum of all values."""
    N = len(barvalues)
    ind = np.arange(N)
    width = 0.35
    plt.title(title)

    # add bars
    sum_values = sum(v for (b, v) in barvalues)
    if sum_values > 0:
        plt.bar(ind, [v * 100 / sum_values for (b, v) in barvalues], width, color='red')
    else:
        plt.bar(ind, [v for (b, v) in barvalues], width, color='red')
    # Final bar is not really visible, increase xlim
    xmin, xmax = plt.xlim()
    plt.xlim(xmin=xmin - 1)
    plt.xlim(xmax=xmax + 1)

    # axis setup
    plt.xticks(ind + width / 2., ["%s\n%s %d" % (b, vlab, v) for (b, v) in barvalues], size='xx-small', rotation=17)
    tick_range = np.arange(0, 110, 10)
    plt.yticks(tick_range, size='xx-small')
    formatter = plt.FixedFormatter([str(x) for x in tick_range])
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().yaxis.grid(which='major')

    savefig_multiple_ext(outdir + "/" + baseout)


def plot_n_save(a, x, y, xlab, ylab, outdir, baseout=False):
    """Plot two columns of a numpy array against each other."""
    # Plot
    plt.plot(a[x], a[y], '.')
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    # Save
    if not baseout:
        baseout = "{x}_vs_{y}".format(x=x, y=y)

    savefig_multiple_ext(outdir + "/" + baseout)


def barplot_read_lcas(a, ranks, title="Barplot read LCA per contig", outdir='.', baseout="barplot_read_lca", read_level="amb", condition=True):
    """Plot the number of times a rank is the LCA of ambiguous reads for a
    contig in a barplot. Expects ranks ordered from lowest rank to highest rank.
    An extra condition over the contig purity numpy array 'a' can be specified
    as an optional argument (e.g. a["max_aln_purity"] < 0.95)."""
    # Compute LCAs
    lcas = []
    has_nolca_yet = True
    # Expect ranks ordered from lowest to highest
    for r in ranks:
        # Taxa should contain all reads if it is LCA
        has_allreads = (a[read_level + '_tot_nr_reads'] == a[read_level + '_nr_reads_' + r]) & ~a.mask[read_level + '_tot_nr_reads']
        is_lca = condition & has_allreads & has_nolca_yet
        lcas.append((r, np.count_nonzero(is_lca)))

        # Determine which contigs have no ambiguous read LCA yet
        has_nolca_yet = has_nolca_yet & ~is_lca

    # Plot the times a rank is LCA per contig
    percent_barplot(lcas, outdir=outdir, baseout=baseout, title=title, vlab="contigs:")


def barplot_read_lcas_per_read(a, ranks, title="Barplot read LCA per read", outdir='.', baseout="barplot_read_lca_per_read", read_level="amb", condition=True):
    lcas = []

    # Expect ranks ordered from lowest to highest
    lower_rank_sum_reads = 0
    for r in ranks:
        num_reads_r = sum(a[~a.mask[read_level + '_tot_nr_reads'] & condition][read_level + '_nr_reads_' + r])
        lcas.append((r, num_reads_r - lower_rank_sum_reads))
        lower_rank_sum_reads = num_reads_r

    # Plot the times a rank is LCA per read
    percent_barplot(lcas, outdir=outdir, baseout=baseout, title=title, vlab="reads:")


def main(genome_contig_cov_tsv, contig_purity_tsv, masmdir):
    # Make plots genome contig coverage
    outdir = masmdir + "/plots"
    make_dir(outdir)
    a = np.genfromtxt(genome_contig_cov_tsv, names=True, dtype=None,
        missing_values=MISSING_VALUE, usemask=True, delimiter=DELIMITER)
    make_dir(outdir)
    plot_n_save(a, "GC_content", "read_cov_ratio", "GC_content of genome", "Read coverage of genome", outdir)
    plot_n_save(a, "GC_content", "genome_contig_cov_ratio", "GC_content of"
            " genome", "Contig coverage ratio of genome", outdir)
    plot_n_save(a, "read_cov_ratio", "genome_contig_cov_ratio", "Read coverage"
            " ratio of genome", "Contig coverage ratio of genome", outdir)

    # Make plots contig purity
    a = np.genfromtxt(contig_purity_tsv, names=True, delimiter=DELIMITER,
                         dtype=None, missing_values=MISSING_VALUE, usemask=True)
    plot_n_save(a, "contig_length", "unamb_read_level_purity", "Contig length",
            "Unambiguous read level purity of contig", outdir)
    plot_n_save(a, "contig_length", "amb_read_level_purity", "Contig length",
            "Ambiguous read level purity of contig", outdir)
    plot_n_save(a, "contig_length", "max_aln_purity", "Contig length",
            "Alignment purity of contig", outdir)
    plot_n_save(a, "unamb_read_level_purity", "max_aln_purity", "Unambiguous"
            " read level purity of contig", "Alignment purity of contig",
            outdir)
    plot_n_save(a, "amb_read_level_purity", "max_aln_purity", "Ambiguous read"
            " level purity of contig", "Alignment purity of contig", outdir)

    is_chimer = a['max_aln_purity'] < 0.95
    is_no_chimer = a['max_aln_purity'] > 0.95
    #ranks_sub = ['strain','species','genus','class','phylum','superkingdom','life']
    ranks_all = ["strain", "species", "genus", "family", "order", "class", "phylum",
            "superphylum", "superkingdom", "life"]
    # Ambiguous reads chimer
#   barplot_read_lcas(a, ranks_sub, title="Ambiguous read LCA for contigs with"
#           " max_aln_purity < 0.95", read_level="amb", condition=is_chimer,
#           outdir=outdir,
#           baseout="barplot_amb_read_lca_subset_ranks_max_aln_purity_lt_95")
    barplot_read_lcas(a, ranks_all, title="Ambiguous read LCA for contigs with"
            " max_aln_purity < 0.95", read_level="amb", condition=is_chimer,
            outdir=outdir,
            baseout="barplot_amb_read_lca_all_ranks_max_aln_purity_lt_95")
    # Unambiguous reads chimer
#   barplot_read_lcas(a, ranks_sub, title="Unambiguous read LCA for contigs with"
#           " max_aln_purity < 0.95", read_level="unamb", condition=is_chimer,
#           outdir=outdir,
#           baseout="barplot_unamb_read_lca_subset_ranks_max_aln_purity_lt_95")
    barplot_read_lcas(a, ranks_all, title="Unambiguous read LCA for contigs with"
            " max_aln_purity < 0.95", read_level="unamb", condition=is_chimer,
            outdir=outdir,
            baseout="barplot_unamb_read_lca_all_ranks_max_aln_purity_lt_95")
    # Amb No chimers
    barplot_read_lcas(a, ranks_all, title="Ambiguous read LCA for contigs with"
            " max_aln_purity > 0.95", read_level="amb", condition=is_no_chimer,
            outdir=outdir,
            baseout="barplot_amb_read_lca_all_ranks_max_aln_purity_gt_95")
    # Unamb no chimer
    barplot_read_lcas(a, ranks_all, title="Unambiguous read LCA for contigs with"
            " max_aln_purity > 0.95", read_level="unamb", condition=is_no_chimer,
            outdir=outdir,
            baseout="barplot_unamb_read_lca_all_ranks_max_aln_purity_gt_95")

    # Amb chimer
    barplot_read_lcas_per_read(a, ranks_all, title="Ambiguous read LCA for contigs with"
            " max_aln_purity < 0.95 per read", read_level="amb", condition=is_chimer,
            outdir=outdir,
            baseout="barplot_amb_read_lca_per_read_all_ranks_max_aln_purity_lt_95")
    # Amb nochimer
    barplot_read_lcas_per_read(a, ranks_all, title="Ambiguous read LCA for contigs with"
            " max_aln_purity > 0.95 per read", read_level="amb", condition=is_no_chimer,
            outdir=outdir,
            baseout="barplot_amb_read_lca_per_read_all_ranks_max_aln_purity_gt_95")
    # Unamb chimer
    barplot_read_lcas_per_read(a, ranks_all, title="Unambiguous read LCA for contigs with"
            " max_aln_purity < 0.95 per read", read_level="unamb", condition=is_chimer,
            outdir=outdir,
            baseout="barplot_unamb_read_lca_per_read_all_ranks_max_aln_purity_lt_95")
    # Unamb no chimer
    barplot_read_lcas_per_read(a, ranks_all, title="Unambiguous read LCA for contigs with"
            " max_aln_purity > 0.95 per read", read_level="unamb", condition=is_no_chimer,
            outdir=outdir,
            baseout="barplot_unamb_read_lca_per_read_all_ranks_max_aln_purity_gt_95")

    # Output plots in HTML
    sdir = os.path.dirname(os.path.realpath(__file__))
    template = open(sdir + '/template/validate-template.html').read()
    with open(masmdir + '/index.html', 'w') as fh:
        fh.write(template.format(asmtsv='asm-stats.tsv',
                                 gcctsv='genome-contig-coverage.tsv',
                                 cptsv='contig-purity.tsv',
                                 plotgc1='plots/GC_content_vs_read_cov_ratio.png',
                                 plotgc2='plots/GC_content_vs_genome_contig_cov_ratio.png',
                                 plotgc3='plots/read_cov_ratio_vs_genome_contig_cov_ratio.png',
                                 plotcp1='plots/contig_length_vs_unamb_read_level_purity.png',
                                 plotcp2='plots/contig_length_vs_amb_read_level_purity.png',
                                 plotcp3='plots/contig_length_vs_max_aln_purity.png',
                                 plotcp4='plots/unamb_read_level_purity_vs_max_aln_purity.png',
                                 plotcp5='plots/amb_read_level_purity_vs_max_aln_purity.png',
                                 plotlcaamb1='plots/barplot_amb_read_lca_all_ranks_max_aln_purity_lt_95.png',
                                 plotlcaamb2='plots/barplot_amb_read_lca_all_ranks_max_aln_purity_gt_95.png',
                                 plotlcaamb3='plots/barplot_amb_read_lca_per_read_all_ranks_max_aln_purity_lt_95.png',
                                 plotlcaamb4='plots/barplot_amb_read_lca_per_read_all_ranks_max_aln_purity_gt_95.png',
                                 plotlcaunamb1='plots/barplot_unamb_read_lca_all_ranks_max_aln_purity_lt_95.png',
                                 plotlcaunamb2='plots/barplot_unamb_read_lca_all_ranks_max_aln_purity_gt_95.png',
                                 plotlcaunamb3='plots/barplot_unamb_read_lca_per_read_all_ranks_max_aln_purity_lt_95.png',
                                 plotlcaunamb4='plots/barplot_unamb_read_lca_per_read_all_ranks_max_aln_purity_gt_95.png',
                                 ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("masmdir", help="Directory of stats output\n")
    args = parser.parse_args()
    masmdir = args.masmdir.rstrip('/')
    main(masmdir + "/genome-contig-coverage.tsv", masmdir + "/contig-purity.tsv", masmdir)
