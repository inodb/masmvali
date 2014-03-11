"""
MASMVALI - Metagenomic ASseMbly VALIdator

Evaluates the quality of a metagenomic assembly given a reference metagenome.
"""
from refgenome import ReferenceSet
import nucmer
import assembly
import argparse
from utils import print_dict2tsv, make_dir, read_2col_table


class AssemblyValidation():
    def __init__(self, *args, **kwargs):
        exp_args = ["bamref", "bamasm", "referenceset", "refstatsfile",
                "refphylfile", "covbedasm", "cut_off", "missing_value",
                "contigs_to_refs_table"]

        for k in kwargs:
            if k not in exp_args:
                raise(Exception("Unexpected keyword received %s" % k))

        self.bamref = kwargs.get("bamref", None)
        self.bamasm = kwargs.get("bamasm", None)
        self.refs = kwargs.get("referenceset", None)
        refstatsfile = kwargs.get("refstatsfile", None)
        refphylfile = kwargs.get("refphylfile", None)
        self.covbedasm = kwargs.get("covbedasm", None)
        self.cut_off = kwargs.get("cut_off", 100)
        self.missing_value = kwargs.get("missing_value", "N/A")
        contigs_to_refs_table = kwargs.get("contigs_to_refs_table", None)
        if contigs_to_refs_table:
            self.contigs_to_refs_dict = read_2col_table(contigs_to_refs_table, sep="\t")
        else:
            self.contigs_to_refs_table = None

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
            self.coords.calc_genome_contig_cov_in_bases(self.refs, cut_off=self.cut_off, contigs_to_refs_dict=self.contigs_to_refs_dict)

        # Print assembly stats
        print_dict2tsv(d=dict(name=self.missing_value, global_purity=self.coords.glob_pur,
                              aln_purity=self.coords.aln_pur, aln_ratio=self.coords.aln_ratio,
                              l50=self.contigs.l50,
                              n50=self.contigs.n50, trim_n=self.contigs.trim_n,
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

    def write_contig_metagenome_cov(self, filename, dec_pow_min=2, dec_pow_max=10, div_per_pow=4):
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
            fh.write("cut_off\tmetagenome_cov_bases\tmetagenome_cov_percentage\n")
            for i in xrange(len(lengthsa)):
                self.coords.calc_genome_contig_cov_in_bases(self.refs,
                        cut_off=int(lengthsa[i]), count_contig_once=True,
                        contigs_to_refs_dict=self.contigs_to_refs_dict)
                fh.write("%d\t%d\t%d\n" % (lengthsa[i], self.coords.ref_sum_cov, (self.coords.ref_sum_cov * 100) / self.refs.mg_length))


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
    val.write_contig_metagenome_cov(masmdir + "/contig-metagenome-coverage.tsv")
