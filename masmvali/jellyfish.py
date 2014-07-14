import sh
import glob
from Bio import SeqIO
from utils import make_dir
from common.cache import property_cached
import os.path


class Jellyfish(object):
    def __init__(self, jf):
        self.jf = jf

    def count_kmer(self, seq):
        awk_count = sh.awk(sh.jellyfish("query", self.jf, "-s", "/dev/fd/0", _in=">seq\n%s" % seq), "{rv+=$2} END {print rv}")
        return int(awk_count.stdout)


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
            sh.jellyfish("count", "-m", kmer, "-s", "10M", "-o",
                    fn.replace(".fa", "-{0}.jf".format(kmer)), fn)

    def get_refs_with_kmer_count_of_seq(self, seq, kmer_size):
        """Check if kmers of given sequence exist in the references"""
        seq = seq.upper()
        refs = {}

        for r in self.records:
            fn = self.records[r].replace(".fa", "-{0}.jf").format(kmer_size)
            if not os.path.exists(fn):
                self._run_kmer(kmer_size)
            refs[r] = Jellyfish(fn).count_kmer(seq)

        return refs
