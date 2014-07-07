"""
MASMVALI - Metagenomic ASseMbly VALIdator

Evaluates the quality of a metagenomic assembly given a reference metagenome.
"""
import os

import numpy as np

from refgenome import ReferenceSet
import nucmer
import assembly
import argparse
from utils import print_dict2tsv, make_dir, read_2col_table

from common.cache import property_cached


class AssemblyValidationParser():
    def __init__(self, directory, missing_value="N/A", sep="\t"):
        self.directory = directory.rstrip('/')
        self.missing_value = missing_value
        self.sep = sep
        self.check_validation_folder()

    def check_validation_folder(self):
        for fp in [self.directory + tsv for tsv in
                ['/asm-stats.tsv', '/genome-contig-coverage.tsv',
                    '/contig-purity.tsv', '/contig-lengths.tsv',
                    '/contig-metagenome-coverage.tsv']]:
            if not os.path.isfile(fp):
                raise(Exception("File {0} not found".format(fp)))

    def _get_numpy_array(self, filepath):
        return np.genfromtxt(filepath, names=True, dtype=None,
                missing_values=self.missing_value, usemask=True, delimiter=self.sep)

    @property_cached
    def asm_stats(self):
        return self._get_numpy_array(self.directory + '/asm-stats.tsv')

    @property_cached
    def genome_contig_coverage(self):
        return self._get_numpy_array(self.directory + '/genome-contig-coverage.tsv')

    @property_cached
    def contig_purity(self):
        return self._get_numpy_array(self.directory + '/contig-purity.tsv')

    @property_cached
    def contig_lengths(self):
        return self._get_numpy_array(self.directory + '/contig-lengths.tsv')

    @property_cached
    def contig_metagenome_coverage_purest_single_aln(self):
        return self._get_numpy_array(self.directory + '/contig-metagenome-coverage-purest-single-aln.tsv')

    @property_cached
    def contig_metagenome_coverage_purest_all_aln(self):
        return self._get_numpy_array(self.directory + '/contig-metagenome-coverage-purest-all-aln.tsv')

    @property_cached
    def contig_metagenome_coverage_all_aln(self):
        return self._get_numpy_array(self.directory + '/contig-metagenome-coverage-all-aln.tsv')

    def get_nx_stats(self, x, cut_off=0):
        cl = self.contig_lengths
        if cl["length"][0] < cut_off:
            raise(Exception("Longest contig shorter than {0}".format(cut_off)))
        # percentage x of sum of contigs longer than cut_off
        sum_cut_off = np.sum(cl["length"][cl["length"] >= cut_off])
        sum_x_cut_off = sum_cut_off * (x / 100.0)
        nx = np.where(cl["sum"] >= sum_x_cut_off)[0][0] + 1
        lx = cl["length"][nx - 1]
        return (lx, nx, sum_x_cut_off, sum_cut_off)

    def get_nmgx_stats(self, x, mg_length):
        """NMGX and LMGX stats as a function of X and MG. To normalize L50 and
        N50 values, one often uses NG50 and NG50 length (LG50) where contigs
        adding up to 50%% of the genome are larger or equal to NG50 length. In
        this case 50 is replaced by X. The function returns (1) NMGX, which is
        the number of contigs adding up to X percent of the length of the
        metagenome (MG). And (2) the shortest length of those contigs i.e.
        LMGX"""
        cl = self.contig_lengths
        if cl["sum"][-1] < (x / 100) * mg_length:
            raise(Exception("Sum of contigs is smaller than x percent of mg_length"))
        nmgx = np.where(cl["sum"] >= mg_length * (x / 100.0))[0][0] + 1
        lmgx = cl["length"][nmgx - 1]
        return (lmgx, nmgx)


class AssemblyValidation():
    def __init__(self, *args, **kwargs):
        self.bamref = kwargs.pop("bamref", None)
        self.bamasm = kwargs.pop("bamasm", None)
        self.refs = kwargs.pop("referenceset", None)
        refstatsfile = kwargs.pop("refstatsfile", None)
        refphylfile = kwargs.pop("refphylfile", None)
        self.covbedasm = kwargs.pop("covbedasm", None)
        self.cut_off = kwargs.pop("cut_off", 100)
        self.missing_value = kwargs.pop("missing_value", "N/A")
        contigs_to_refs_table = kwargs.pop("contigs_to_refs_table", None)
        if contigs_to_refs_table:
            self.contigs_to_refs_dict = read_2col_table(contigs_to_refs_table, sep="\t")
        else:
            self.contigs_to_refs_table = None
        if len(kwargs) > 0:
            raise(Exception("Unexpected keyword argument found %s" % kwargs))

        self.asmfa = args[0]
        self.coordsfile = args[1]

        # Determine reference metagenome
        if not self.refs:
            if refphylfile and refstatsfile:
                self.refs = ReferenceSet(refphylfile, refstatsfile, contigs_to_refs_dict=self.contigs_to_refs_dict)
            else:
                raise(Exception("Either supply refphylfile and refstatsfile or an existing reference set"))

        # Determine read-contig and read-reference mappings from bam files if
        # available
        if self.bamref and self.bamasm and self.covbedasm:
            self.reads, self.contigs = assembly.get_read_contig_mappings(self.bamref,
                                                                         self.bamasm,
                                                                         self.refs,
                                                                         self.asmfa,
                                                                         self.cut_off)
            self.contigs.parse_cov_bed(self.covbedasm)
        else:
            self.contigs = assembly.ContigDict()
            self.contigs.parse_fasta_lengths(self.asmfa)

        # Calculate stats from nucmer alignments
        self.coords = nucmer.Coords(self.coordsfile)
        self.coords.calc_max_aln_purity_per_contig(self.contigs, self.cut_off)

    def write_contig_purity(self, filename):
        with open(filename, "w") as fh:
            fh.write(self.contigs.itervalues().next().str_stats_header() + "\n")
            for c in self.contigs.itervalues():
                fh.write(c.str_stats() + "\n")

    def write_general_stats(self, filename):
        if not hasattr(self.coords, "q_aln_bases"):
            self.coords.calc_alignedbases_per_contig(self.contigs)

        # Get all reference lengths
        if not hasattr(self.coords, "ref_sum_bases"):
            self.coords.calc_genome_contig_cov_in_bases(self.refs,
                    cut_off=self.cut_off, count_contig_once=True,
                    only_purest_alignments=True,
                    contigs_to_refs_dict=self.contigs_to_refs_dict)

        # Print assembly stats
        nmg50, lmg50 = self.contigs.get_ng50_and_lg50(self.refs.mg_length) or (self.missing_value, self.missing_value)
        print_dict2tsv(d=dict(name=self.missing_value, global_purity=self.coords.glob_pur,
                              aln_purity=self.coords.aln_pur, aln_ratio=self.coords.aln_ratio,
                              l50=self.contigs.l50,
                              n50=self.contigs.n50, nmg50=nmg50, lmg50=lmg50, trim_n=self.contigs.trim_n,
                              max_contig_length=self.contigs.max_length,
                              cut_off=self.cut_off,
                              trim_n_mapping=sum([hasattr(c, "aln_bases") for c in self.contigs.itervalues()]),
                              sum_ref_lengths=self.refs.mg_length,
                              sum_purest_bases=self.coords.sum_purest_bases,
                              metagenome_cov=float(self.coords.ref_sum_cov) / self.refs.mg_length,
                              sum_bases=self.contigs.totbases, asm_type=self.missing_value,
                              kmer_type=self.missing_value, kmer_size=self.missing_value,
                              kmin=self.missing_value,
                              kmax=self.missing_value), filepath=filename)

    def write_genome_contig_cov(self, filename):
        if not hasattr(next(iter(self.refs.refs.values())), "contig_cov"):
            self.coords.calc_genome_contig_cov_in_bases(self.refs, cut_off=self.cut_off, contigs_to_refs_dict=self.contigs_to_refs_dict)
        with open(filename, "w") as fh:
            fh.write("genome\tgenome_contig_cov_bases\tgenome_length\tgenome_contig_cov_ratio\tGC_content\tread_cov_ratio\tread_cov_mean\n")
            for r in set(self.refs.refs.values()):
                fh.write("%s\t%i\t%i\t%f\t%f\t%f\t%f\n" %
                        (r.phyl["strain"],
                         getattr(r, "contig_cov", 0), r.length, getattr(r, "contig_cov",
                         0) / float(r.length), r.gc_content, r.ratio_covered, r.cov))

    def write_contig_lengths(self, filename):
        sum_lens = 0
        with open(filename, "w") as fh:
            fh.write("name\tlength\tsum\n")
            for name in sorted(self.contigs, key=lambda x: self.contigs.get(x).length, reverse=True):
                sum_lens += self.contigs[name].length
                fh.write("%s\t%d\t%d\n" % (name, self.contigs[name].length, sum_lens))

    def write_contig_metagenome_cov(self, filename, dec_pow_min=2, dec_pow_max=10, div_per_pow=4, incl_alignments='purest-single'):
        """Write metagenome coverage to a file. The different alignments that
        can be used are 'purest-single' for one purest alignment per contig,
        'purest-all' for all purest alignments per contig 'all' for all
        alignments regardless of purity."""
        if incl_alignments == 'purest-single':
            count_contig_once = True
            only_purest_alignments = True
        elif incl_alignments == 'purest-all':
            count_contig_once = False
            only_purest_alignments = True
        elif incl_alignments == 'all':
            count_contig_once = False
            only_purest_alignments = False
        else:
            raise(Exception("Invalid value for incl_alignments, should be purest-single, purest-all or all: " + incl_alignments))

        # Only a subset of the contig lengths are used to determine metagenome
        # coverage, because using all takes too much time. Make list of
        # specified powers of 10, divide by the number of divisions one wants
        # to have for each power. Sort and reverse.
        lengthsa = \
            [item for sublist in map(lambda x: (range(10 ** x, 10 ** (x + 1), (10 ** (x + 1) - 10 ** x) / div_per_pow)), range(dec_pow_min, dec_pow_max)) for item in sublist] + \
            [10 ** dec_pow_max] + [self.contigs.l50]
        lengthsa.sort()
        lengthsa.reverse()
        with open(filename, "w") as fh:
            fh.write("cut_off\tmetagenome_cov_bases\tmetagenome_cov_percentage\tassembly_bases\n")
            for i in xrange(len(lengthsa)):
                self.coords.calc_genome_contig_cov_in_bases(self.refs,
                    cut_off=int(lengthsa[i]), count_contig_once=count_contig_once,
                    only_purest_alignments=only_purest_alignments,
                    contigs_to_refs_dict=self.contigs_to_refs_dict)
                fh.write("%d\t%d\t%d\t%d\n" % (lengthsa[i],
                    self.coords.ref_sum_cov, (self.coords.ref_sum_cov * 100) /
                    self.refs.mg_length,
                    sum([c.length for c in self.contigs.itervalues() if c.length >= lengthsa[i]])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("coords", help="Output from nucmer show-coords")
    parser.add_argument("refstats", help="Reference stats tsv file\n")
    parser.add_argument("refphyl", help="Reference Phylogeny tsv file\n")
    parser.add_argument(
        "contigs", help="Contigs of the assembly in fasta format\n")
    parser.add_argument("masmdir", help="Directory to output stats\n")
    parser.add_argument(
        "--bamref", default=None, help="BAM file of the reads mapped against the reference\n")
    parser.add_argument(
        "--bamasm", default=None, help="BAM file of the reads mapped against the contigs\n")
    parser.add_argument(
        "--covbedasm", default=None, help="coverageBed output of BAM file of the reads mapped against the contigs\n")
    parser.add_argument(
        "--contigs_to_refs_table", default=None, help="Table in tsv format,"
        " first column contig names, second column reference names, no header\n")
    args = parser.parse_args()

    if (args.bamref or args.bamasm or args.covbedasm) and not (args.bamref and args.bamasm and args.covbedasm):
        raise(Exception("All three options --bamref, --bamasm and --covbedasm required if read purity"
                        " computation is wanted"))

    # Calculate stats
    val = AssemblyValidation(args.contigs, args.coords, refphylfile=args.refphyl,
                             refstatsfile=args.refstats, bamref=args.bamref,
                             bamasm=args.bamasm, covbedasm=args.covbedasm, contigs_to_refs_table=args.contigs_to_refs_table)

    # Print output
    masmdir = args.masmdir.rstrip('/')
    make_dir(masmdir)
    val.write_general_stats(masmdir + "/asm-stats.tsv")
    val.write_genome_contig_cov(masmdir + "/genome-contig-coverage.tsv")
    val.write_contig_purity(masmdir + "/contig-purity.tsv")
    val.write_contig_lengths(masmdir + "/contig-lengths.tsv")
    #val.write_contig_metagenome_cov(masmdir + "/contig-metagenome-coverage-purest-single-aln.tsv", incl_alignments='purest-single')
    #val.write_contig_metagenome_cov(masmdir + "/contig-metagenome-coverage-purest-all-aln.tsv", incl_alignments='purest-all')
    val.write_contig_metagenome_cov(masmdir + "/contig-metagenome-coverage-all-aln.tsv", incl_alignments='all')
