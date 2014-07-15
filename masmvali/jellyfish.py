from __future__ import print_function
import sh
import glob
from Bio import SeqIO
from utils import make_dir
from common.cache import property_cached
import os.path
import sys
from subprocess import Popen, PIPE


class Jellyfish(object):
    def __init__(self, jf):
        self.jf = jf

    def count_existent_kmers(self, seq):
        """Checks whether the kmers in the given sequence exist in the
        jellyfish database. The number of kmers that exist are returned (could
        be multiple times the same kmer)."""
        #awk_count = sh.awk(sh.jellyfish("query", self.jf, "-s", "/dev/fd/0", _in=">seq\n%s\n" % seq, _piped=True), "{rv+=$2} END {print rv}")
        #awk_count = sh.jellyfish("query", self.jf, "-s", "/dev/fd/0", '|', "{rv+=$2} END {print rv}", _in=">seq\n%s\n" % seq)
        p1 = Popen(["jellyfish", "query", self.jf, "-s", "/dev/fd/0"], stdout=PIPE, stdin=PIPE)
        p2 = Popen(["awk", "BEGIN {rv=0} { if ($2>0) {rv+=1} } END {print rv}"], stdin=p1.stdout, stdout=PIPE)
        #p1 = Popen(["jellyfish", "query", self.jf, "-s", "/dev/fd/0"], stdout=PIPE, stdin=PIPE)
        p1.communicate(input=">seq\n%s\n" % seq)
        p2rv = p2.communicate()
        awk_count = p2rv[0]
        return int(awk_count)

    def count_non_existent_kmers(self, seq):
        """Checks whether the kmers in the given sequence exist in the
        jellyfish database. The number of kmers that don't exist are returned
        (could be multiple times the same kmer)."""
        p1 = Popen(["jellyfish", "query", self.jf, "-s", "/dev/fd/0"], stdout=PIPE, stdin=PIPE)
        p2 = Popen(["awk", "BEGIN {rv=0} { if ($2==0) {rv+=1} } END {print rv}"], stdin=p1.stdout, stdout=PIPE)
        p1.communicate(input=">seq\n%s\n" % seq)
        p2rv = p2.communicate()
        awk_count = p2rv[0]
        return int(awk_count)


class JellyfishRef(object):
    def __init__(self, reffasta, folder):
        self.jf_folder = folder.rstrip("/")
        self.reffasta = reffasta
        make_dir(self.jf_folder)
        self._init_ref(reffasta)

    def _init_ref(self, reffasta):
        """Split the reference fasta over multiple fasta files for querying
        with Jellyfish."""
        self.records = {}

        # split the fasta files
        for record in SeqIO.parse(reffasta, "fasta"):
            filename = self.jf_folder + "/" + record.id.replace(" ",
                    "_").replace("/", "_") + ".fa"
            with open(filename, "w") as fh:
                fh.write(">%s\n%s" % (record.id, record.seq))
            self.records[record.id] = filename

    def _run_kmer(self, kmer):
        """Run kmer counter for given kmer size on all references"""
        for fn in self.records.itervalues():
            sh.jellyfish("count", "-C", "-m", kmer, "-s", "10M", "-o",
                    fn.replace(".fa", "-{0}.jf".format(kmer)), fn)
        # one for all genomes combined (count both strands -C)
        sh.jellyfish("count", "-C", "-m", kmer, "-s", "10M", "-o",
                "{0}/ref-{1}.jf".format(self.jf_folder, kmer), *self.records.itervalues())

    def get_refs_with_kmer_count_of_seq(self, seq, kmer_size):
        """Check if kmers of given sequence exist in the references"""
        seq = seq.upper()
        refs = {}

        for r in self.records:
            fn = self.records[r].replace(".fa", "-{0}.jf").format(kmer_size)
            if not os.path.exists(fn):
                self._run_kmer(kmer_size)
            refs[r] = Jellyfish(fn).count_existent_kmers(seq)

        refs["invalid"] = Jellyfish("{0}/ref-{1}.jf".format(self.jf_folder,
            kmer_size)).count_non_existent_kmers(seq)

        return refs

    def write_refs_kmer_count_of_fasta(self, fafile, kmer_size, output_file):
        with open(output_file, "w") as fh:
            fh.write("name\tlength\tinvalid\t{0}\n".format("\t".join([r for r in sorted(self.records)])))

            for record in SeqIO.parse(fafile, "fasta"):
                kmer_counts = self.get_refs_with_kmer_count_of_seq(record.seq, kmer_size)
                fh.write("{0}\t{1}\t{2}\t{3}\n".format(record.id, len(record.seq), kmer_counts["invalid"], "\t".join([str(kmer_counts[r]) for r in sorted(self.records)])))
