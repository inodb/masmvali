from refgenome import ReferenceSet
import nucmer
import assembly
import argparse
from utils import print_dict2tsv, make_dir


class AssemblyValidation():
    def __init__(self, refphylfile, refstatsfile, asmfa, nucmercoords,
                 bamref=None, bamasm=None, covbedasm=None, missing_value="N/A", cut_off=100):
        self.bamref = bamref
        self.bamasm = bamasm
        self.covbedasm = covbedasm
        self.asmfa = asmfa
        self.cut_off = cut_off
        self.missing_value = missing_value

        # Determine reference metagenome
        self.refs = ReferenceSet(refphylfile, refstatsfile)

        # Determine read-contig and read-reference mappings from bam files if
        # available
        if bamref and bamasm and covbedasm:
            self.reads, self.contigs = assembly.get_read_contig_mappings(bamref,
                bamasm, self.refs, asmfa, self.cut_off)
            self.contigs.parse_cov_bed(covbedasm)
        else:
            self.contigs = assembly.ContigDict()
            self.contigs.parse_fasta_lengths(asmfa)

        # Calculate stats from nucmer alignments
        self.coords = nucmer.Coords(nucmercoords)
        self.coords.calc_max_aln_purity_per_contig(self.contigs, self.cut_off)

    def write_contig_purity(self, filename, sep="\t"):
        with open(filename, "w") as fh:
            fh.write(self.contigs.itervalues().next().str_stats_header() + "\n")
            for c in self.contigs.itervalues():
                fh.write(c.str_stats() + "\n")

    def write_general_stats(self, filename, sep="\t"):
        q_aln_bases = self.coords.calc_alignedbases_per_contig()
        sum_purest_bases = sum(
            getattr(c, "max_aln_purity", 0) * c.length for c in
            self.contigs.itervalues())
        sum_aln_bases = sum(q_aln_bases.itervalues())

        # Alignment purity, how pure are all the alignments
        aln_pur = float(sum_purest_bases) / sum_aln_bases
        # Global purity
        glob_pur = float(sum_purest_bases) / self.contigs.totbases
        aln_ratio = float(sum_aln_bases) / self.contigs.totbases

        # Get all reference lengths
        ref_sum_bases = sum(r.length for r in set(self.refs.refs.values()))
        if not hasattr(next(iter(self.refs.refs.values())), "contig_cov"):
            self.coords.calc_genome_contig_cov_in_bases(self.refs.refs)
        ref_sum_cov = sum([getattr(r, "contig_cov", 0) for r in
                            set(self.refs.refs.values())])

        # Print assembly stats
        print_dict2tsv(d=dict(name=self.missing_value, global_purity=glob_pur,
                            aln_purity=aln_pur, aln_ratio=aln_ratio,
                            l50=self.contigs.l50,
                            n50=self.contigs.n50, trim_n=self.contigs.trim_n,
                            max_contig_length=self.contigs.max_length,
                            cut_off=self.cut_off, trim_n_mapping=len(q_aln_bases),
                            sum_ref_lengths=ref_sum_bases,
                            sum_purest_bases=sum_purest_bases,
                            metagenome_cov=float(ref_sum_cov) / ref_sum_bases,
                            sum_bases=self.contigs.totbases, asm_type=self.missing_value,
                            kmer_type=self.missing_value, kmer_size=self.missing_value,
                            kmin=self.missing_value,
                            kmax=self.missing_value), filepath=filename)

    def write_genome_contig_cov(self, filename, sep="\t"):
        if not hasattr(next(iter(self.refs.refs.values())), "contig_cov"):
            self.coords.calc_genome_contig_cov_in_bases(self.refs.refs)
        with open(filename, "w") as fh:
            fh.write("genome\tgenome_contig_cov_bases\tgenome_length\tgenome_contig_cov_ratio\tGC_content\tread_cov_ratio\tread_cov_mean\n")
            for r in set(self.refs.refs.values()):
                fh.write("%s\t%i\t%i\t%f\t%f\t%f\t%f\n" % (r.phyl["strain"],
                    getattr(r, "contig_cov", 0), r.length, getattr(r, "contig_cov",
                    0) / float(r.length), r.gc_content, r.ratio_covered, r.cov))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
    args = parser.parse_args()

    if (args.bamref or args.bamasm or args.covbedasm) and not (args.bamref and args.bamasm and args.covbedasm):
        raise(Exception("All three options --bamref, --bamasm and --covbedasm required if read purity"
                        " computation is wanted"))

    # Calculate stats
    val = AssemblyValidation(args.refphyl, args.refstats, args.contigs,
    args.coords, args.bamref, args.bamasm, args.covbedasm)

    # Print output
    masmdir = args.masmdir.rstrip('/')
    make_dir(masmdir)
    val.write_general_stats(masmdir + "/asm-stats.tsv")
    val.write_genome_contig_cov(masmdir + "/genome-contig-coverage.tsv")
    val.write_contig_purity(masmdir + "/contig-purity.tsv")
