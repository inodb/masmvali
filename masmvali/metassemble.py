import warnings
import re
import os
import numpy as np

from refgenome import ReferenceSet
from validation import AssemblyValidation, AssemblyValidationParser
from utils import read_2col_table


class MetAssemble():
    """Class to handle results from MetAssemble pipeline."""
    # strategies without kmers
    all_strategies_alpha = \
        ['velvetnoscaf',
         'velvetscaf',
         'velvetnoscafminimus2',
         'velvetnoscafnewbler',
         'velvetnoscafbambus2',
         'velvetnoscafminimus2bambus2',
         'velvetnoscafnewblerbambus2',
         'metavelvetnoscaf',
         'metavelvetscaf',
         'metavelvetnoscafminimus2',
         'metavelvetnoscafnewbler',
         'metavelvetnoscafbambus2',
         'metavelvetnoscafminimus2bambus2',
         'metavelvetnoscafnewblerbambus2',
         'raynoscaf',
         'rayscaf',
         'raynoscafminimus2',
         'raynoscafnewbler',
         'raynoscafbambus2',
         'raynoscafminimus2bambus2',
         'raynoscafnewblerbambus2']

    fasta_postfixes = \
        ['/noscaf/noscaf_X/ma-contigs.fa',
         '/scaf/scaf_X/ma-scaffolds.fa',
         '/noscaf/minimus2/ma-merge.fa',
         '/noscaf/newbler/ma-merge.fa',
         '/noscaf/noscaf_X/bambus2/bambus2.scaffold.linear.fasta',
         '/noscaf/minimus2/bambus2/bambus2.scaffold.linear.fasta',
         '/noscaf/newbler/bambus2/bambus2.scaffold.linear.fasta']

    # fasta file paths dictionary
    fasta_fpd = \
        dict(zip(all_strategies_alpha,
                 [item for sublist in
                 [['/velvet' + pf for pf in fasta_postfixes],
                  ['/metavelvet' + pf for pf in fasta_postfixes],
                  ['/ray' + pf for pf in fasta_postfixes]]
                  for item in sublist]))

    def __init__(self, metassembledir, refphylfile, refstatsfile, contigs_to_refs_file, sep="\t"):
        self.metassembledir = metassembledir.rstrip('/')
        self.refphylfile = refphylfile
        self.refstatsfile = refstatsfile
        self.contigs_to_refs_file = contigs_to_refs_file
        self.contigs_to_refs_dict = read_2col_table(contigs_to_refs_file, sep=sep)
        self.refs = ReferenceSet(refphylfile, refstatsfile, contigs_to_refs_dict=self.contigs_to_refs_dict)

    def alphafy_strat(self, strategy):
        alpha_strat = re.sub("scaf[0-9]{2}", "scaf", strategy)
        if alpha_strat not in self.fasta_fpd:
            raise(Exception("Invalid strategy: " + alpha_strat))
        else:
            return alpha_strat

    def has_kmer(strategy):
        # Should consist of 2 numbers and should not be 22 (this is minimus2bambus2)
        return len(re.sub("[a-z]", "", strategy)) > 1 and int(re.sub("[a-z]", "", strategy)) != 22

    def get_kmer(self, strategy):
        # Only take first two numbers, the 2 in bambus2 can add a number
        numbers = re.sub("[a-z]", "", strategy)
        return int(numbers[0:2])

    def get_asm_strategies(self, min_kmer, max_kmer, stepsize):
        add_kmers = lambda name: [name + str(n) for n in range(min_kmer, max_kmer, stepsize)]

        asm_sets = {}
        asm_sets["contig"] = add_kmers("velvetnoscaf") + add_kmers("raynoscaf") + add_kmers("metavelvetnoscaf")
        asm_sets["scaf"] = add_kmers("velvetscaf") + add_kmers("rayscaf") + add_kmers("metavelvetscaf")
        asm_sets['contig-merge'] = [n for n in self.all_strategies_alpha if n.endswith('newbler') or n.endswith('minimus2')]
        asm_sets['contig-scaf'] = [a + "bambus2" for a in asm_sets["contig"]]
        asm_sets['contig-merge-scaf'] = [a + "bambus2" for a in asm_sets["contig-merge"]]

        return asm_sets

    def get_metassemble_assembly_path(self, strategy):
        d = self.metassembledir

        # Take only name without kmer number in it
        fp = self.fasta_fpd[self.alphafy_strat(strategy)]

        # Change X to requested kmer size
        if re.search("X", fp):
            fp = re.sub("X", str(self.get_kmer(strategy)), fp)

        fp = d + '/assemblies' + fp
        if not os.path.isfile(fp):
            warnings.warn("Strategy " + strategy + "'s fasta file does not exist: " + fp)

        return fp

    def get_metassemble_assembly_validation(self, strategy):
        fp = self.get_metassemble_assembly_path(strategy)
        return AssemblyValidation(fp, "/".join(fp.split("/")[:-1]) + '/val/nucmer.coords', referenceset=self.refs)

    def get_metassemble_assembly_validation_parser(self, strategy, missing_value="N/A"):
        fp = self.get_metassemble_assembly_path(strategy)
        return AssemblyValidationParser("/".join(fp.split("/")[:-1]) + '/val', missing_value)

    def get_metassemble_assemstats(self, strategy, MISSING_VALUE="N/A", DELIMITER="\t"):
        """returns numpy array of asm-stats.tsv"""
        fp = self.get_metassemble_assembly_path(strategy)
        a = np.genfromtxt("/".join(fp.split("/")[:-1]) + '/val/asm-stats.tsv',
                          names=True, dtype=None, missing_values=MISSING_VALUE, usemask=True,
                          delimiter=DELIMITER)
        return a
