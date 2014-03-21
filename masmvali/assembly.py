from collections import namedtuple, Counter, defaultdict, MutableMapping
import re

import pysam

from refgenome import EqualNameImpliesEquality, Reference


class Read(EqualNameImpliesEquality):
    def __init__(self, name):
        self.name = name
        self.references = set()
        self.contigs = set()

    def add_ref(self, ref):
        self.references.add(ref)

    def add_contig(self, contig):
        self.contigs.add(contig)

    def num_ref(self):
        return(len(self.references))


class ReadDict():
    ReadStats = namedtuple('ReadStats', 'dom_strain, dom_nr_reads, read_level_purity')

    def __init__(self):
        self.noref_reads = set()
        self.unamb_reads = set()
        self.amb_reads = set()

        self.amb_stats = None
        self.unamb_stats = None

    def __len__(self):
        return(len(self.noref_reads) + len(self.unamb_reads) + len(self.amb_reads))

    def add(self, read):
        num_refs = read.num_ref()

        if num_refs == 0:
            self.noref_reads.add(read)
        elif num_refs == 1:
            self.unamb_reads.add(read)
        else:
            self.amb_reads.add(read)

    def get_unamb_stats(self):
        if self.unamb_stats is not None:
            return(self.unamb_stats)
        else:
            unamb_count = Counter()
            for ur in self.unamb_reads:
                unamb_count.update(ur.references)
            mc = unamb_count.most_common(1)

            if len(mc) == 0:
                return None

            self.unamb_stats = self.ReadStats(dom_strain=mc[0][0], dom_nr_reads=mc[0][1],
                                              read_level_purity=float(mc[0][1])
                                                  / len(self.unamb_reads))
            return(self.unamb_stats)

    def get_amb_stats(self):
        # TODO: fix this weird lazy recomputing stuff to a better solution
        if self.amb_stats is not None:
            return(self.amb_stats)
        else:
            # Determine dominant strain including ambiguous read-ref hit counts
            amb_count = Counter()
            for ur in self.unamb_reads:
                amb_count.update(ur.references)
            for ar in self.amb_reads:
                amb_count.update(ar.references)
            mc = amb_count.most_common(1)

            if len(mc) == 0:
                return None

            self.amb_stats = self.ReadStats(dom_strain=mc[0][0], dom_nr_reads=mc[0][1],
                                     read_level_purity=float(mc[0][1]) /
                                     (len(self.unamb_reads) + len(self.amb_reads)))
            return(self.amb_stats)

    def get_unamb_nr_reads_per_taxon(self):
        #TODO: pretty similar to get_amb_nr_reads, can't they be generalized?
        # Maybe have read become an UnAmbiguous read type and override get_lca,
        # some for get_stats
        tax_count = defaultdict(int)
        if self.get_unamb_stats() is not None:
            for ur in self.unamb_reads:
                tax_count[iter(ur.references).next().get_lca(self.unamb_stats.dom_strain)[0]] += 1

            # Add LCAs of lower taxonomic ranks to counts of higher ones (i.e. make
            # the counts ascending with increasing taxonomic rank):
            tot = 0
            for t in Reference.TAX_LVLS:
                tax_count[t] += tot
                tot = tax_count[t]

            return(tax_count)
        else:
            return(dict([(t, 0) for t in Reference.TAX_LVLS]))

    def get_amb_nr_reads_per_taxon(self):
        amb_stats = self.get_amb_stats()
        if amb_stats is not None:
            tax_count = defaultdict(int)
            for ur in self.unamb_reads:
                tax_count[iter(ur.references).next().get_lca(self.unamb_stats.dom_strain)[0]] += 1
            # Determine tax count based on the dominant strain. Every ambiguous
            # read is assigned to the reference with the highest LCA. TODO: In
            # case the dominant strain is also part of the references of the
            # ambiguous reads, don't we want it to pick the lowest LCA in that
            # case?
            for ar in self.amb_reads:
                tax_count[amb_stats.dom_strain.get_highest_lca(ar.references)[0]] += 1

            # Add LCAs of lower taxonomic ranks to counts of higher ones (i.e. make
            # the counts ascending with increasing taxonomic rank):
            tot = 0
            for t in Reference.TAX_LVLS:
                tax_count[t] += tot
                tot = tax_count[t]

            return(tax_count)
        else:
            return(self.get_unamb_nr_reads_per_taxon())

    @staticmethod
    def get_stats_header():
        rv = []

        for prefix in ("unamb_", "amb_"):
            rv.extend([prefix + col for col in ["dominant_strain",
                "dom_nr_reads", "read_level_purity", "tot_nr_reads"] +
                        ["nr_reads_" + t for t in Reference.TAX_LVLS]])

        return rv

    @staticmethod
    def str_stats_header(sep="\t"):
        return(sep.join(ReadDict.get_stats_header()))

    def str_stats(self, sep="\t", missing_value="N/A"):
        rv = ""

        unamb_stats = self.get_unamb_stats()
        if unamb_stats is not None:
            rv += sep.join([str(s) for s in unamb_stats])
            rv += sep + str(len(self.unamb_reads))

            unamb_taxon_counts = self.get_unamb_nr_reads_per_taxon()
            for t in Reference.TAX_LVLS:
                rv += sep + str(unamb_taxon_counts[t])
        else:
            rv += sep.join([missing_value] * (len(ReadDict.get_stats_header()) / 2))

        amb_stats = self.get_amb_stats()
        if amb_stats is not None:
            for s in amb_stats:
                rv += sep + str(s)
            # TODO: This is a bit weird because the outputted ambiguous reads
            # are the ambiguous reads incl unambiguous reads whereas internal
            # they are called amb_reads and unamb_reads. Should reconsider this
            # design.
            rv += sep + str(max(len(self.amb_reads), len(self.unamb_reads)))

            amb_taxon_counts = self.get_amb_nr_reads_per_taxon()
            for t in Reference.TAX_LVLS:
                rv += sep + str(amb_taxon_counts[t])
        else:
            rv += sep + sep.join([missing_value] * (len(ReadDict.get_stats_header()) / 2))

        return(rv)


class Contig(EqualNameImpliesEquality):
    def __init__(self, name):
        self.name = name

    def add_read(self, read):
        try:
            self.reads.add(read)
        except AttributeError:
            self.reads = ReadDict()
            self.reads.add(read)

    def str_stats_header(self, sep="\t"):
        return(sep.join(["contig", "contig_length", "max_aln_strain", "max_aln_purity", "cov_mean"] + ReadDict.get_stats_header()))

    def str_stats(self, sep="\t", missing_value="N/A"):
        rv = ""

        rv += sep.join([self.name, str(self.length),
            str(getattr(self, "max_aln_strain", missing_value)),
            str(getattr(self, "max_aln_purity", missing_value)),
            str(getattr(self, "cov_mean", missing_value))])
        rv += sep
        if hasattr(self, "reads"):
            rv += self.reads.str_stats(missing_value=missing_value)
        else:
            rv += sep.join([missing_value] * len(ReadDict.get_stats_header()))

        return(rv)


class ContigDict(MutableMapping):
    def __init__(self):
        self.contigs = {}

    def __getitem__(self, key):
        return(self.contigs[key])

    def __setitem__(self, key, value):
        #TODO:
        raise(Exception("ContigDict can only add new contigs through parse_fasta_lengths()"))

    def __delitem__(self, key):
        del self.contigs[key]

    def __len__(self):
        return(len(self.contigs))

    def __iter__(self):
        return(iter(self.contigs))

    def parse_fasta_lengths(self, fafile, cut_off=100):
        self.l50, self.n50, self.totbases, self.max_length, self.trim_n, fa_contig_lengths = asm_stats_fasta(fafile, cut_off)
        for factot in fa_contig_lengths:
            #TODO: temp fix with splitting factot because nucmer doesnt store entire contig name, bam prolly neither
            fac = factot.split()[0]
            try:
                c = self.contigs[fac]
            except KeyError:
                c = Contig(fac)
                self.contigs[fac] = c
            c.length = fa_contig_lengths[factot]

    def get_ng50_and_lg50(self, genome_length):
        sum_lens = 0
        nr = 0
        for name in sorted(self, key=lambda x: self.get(x).length, reverse=True):
            nr += 1
            sum_lens += self.get(name).length
            if sum_lens >= genome_length:
                return (nr, self.get(name).length)
        return None

    def parse_cov_bed(self, covbedasm):
        """Take the default coverageBed output and store coverage mean for each contig in self."""
        cov_means = defaultdict(int)
        # Parse the coverageBed output file, 0 name, 1 cov, 2 nr of bases with
        # that cov, 3 total nr of bases in contig, 4 ratio of total bases in
        # contig
        for line in open(covbedasm):
            cols = line.split()
            cov_means[cols[0]] = cov_means[cols[0]] + int(cols[1]) * float(cols[4])

        # Copy to self
        for cm in cov_means:
            try:
                self.contigs[cm].cov_mean = cov_means[cm]
            except KeyError:
                if cm == "genome":
                    self.cov_mean = cov_means["genome"]
                else:
                    raise(Exception("Contig %s doesn't exist. ContigDict can only add new contigs through parse_fasta_lengths()" % cols[0]))


def get_read_contig_mappings(bamref, bamasm, refs, asmfa, cut_off=100):
    """Takes a bam file of reads mapped to the references and a bam file of
    reads mapped against the assembly and returns a dictionary of reads indexed
    by name and contigs indexed by name.

    bamref -- bamfile of reads mapped to reference
    bamasm -- bamfile of the same reads mapped to contigs/scaffolds of an
              assembly
    refs   -- A ReferenceSet
    asmfa  -- Fasta file of the assembly"""
    reffile = pysam.Samfile(bamref, "rb")
    reads = {}
    for record in reffile:
        if record.mapq > 0:
            try:
                r = reads[record.qname]
            except KeyError:
                r = Read(record.qname)
                reads[record.qname] = r
            ref = refs.get(reffile.getrname(record.tid))
            r.add_ref(ref)

    contigs = ContigDict()
    contigs.parse_fasta_lengths(asmfa)
    asmfile = pysam.Samfile(bamasm, "rb")
    for record in asmfile:
        if record.mapq > 0 and record.qlen >= cut_off:
            # All contigs should be in asmfa
            # TODO: use try/except return sensible error msg
            c = contigs.get(asmfile.getrname(record.tid))
            try:
                r = reads[record.qname]
            except KeyError:
                r = Read(record.qname)
                reads[record.qname] = r
            c.add_read(r)
            r.add_contig(c)

    return reads, contigs


def asm_stats_fasta(fafile, cut_off=100):
    """Return L50 of an assembly and the total number of bases. 50% of all
    bases in the assembly are located in contigs equal or larger than l50."""
    name = None
    all_contig_lengths = {}

    # Determine lengths of all contigs
    for line in open(fafile):
        if line[0] == '>':
            name = line[1:-1]
            all_contig_lengths[name] = 0
        else:
            all_contig_lengths[name] += len("".join(re.findall(
                "[ACGTacgt]+", line)))

    # Determine l50, 50% of the bases in the fasta file are located in contigs => l50
    sorted_contigs = sorted(filter(lambda x: x >= cut_off,
                                   all_contig_lengths.values()),
                            reverse=True)
    cumsum = 0
    totbases = sum(sorted_contigs)
    max_length = sorted_contigs[0]
    for i in range(len(sorted_contigs)):
        cumsum += sorted_contigs[i]
        if cumsum >= totbases / 2:
            l50 = sorted_contigs[i]
            n50 = i + 1
            break

    return(l50, n50, totbases, max_length, len(sorted_contigs), all_contig_lengths)
