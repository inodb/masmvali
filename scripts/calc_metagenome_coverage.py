import sys
import os
import argparse

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../masmvali/')
import nucmer


from Bio import SeqIO


def calc_metagenome_coverage(coords, reffa, cut_off):
    refs = get_gc_and_len_dict(reffa)
    ccursor = nucmer.get_coords_db_cursor(coords)
    gcc = nucmer.calc_genome_contig_cov_in_bases(ccursor, cut_off, count_contig_once=True)
    for k in gcc:
        if k not in refs:
            raise(Exception('Reference in coords file not in reference fasta: %s' % k))
        else:
            refs[k]["contig_cov"] = gcc[k]
    return refs


def print_metagenome_coverage(coords, reffa, cut_off):
    refs = calc_metagenome_coverage(coords, reffa, cut_off)

    # per genome
    print "%s\t%s\t%s\t%s" % ('genome', 'length', 'contig_cov', 'contig_cov_percentage')
    for k in refs:
        print "%s\t%i\t%i\t%.2f" % (k,
            refs[k]["length"],
            refs[k].get("contig_cov", 0),
            refs[k]["contig_cov"] * 100.0 / refs[k]["length"] if "contig_cov" in refs[k] else 0)
    
    # entire metagenome
    ref_sum_bases = sum([refs[k]["length"] for k in refs])
    contig_cov_sum_bases = sum([refs[k].get("contig_cov", 0) for k in refs])
    print "%s\t%i\t%i\t%.2f" % ('metagenome',
        ref_sum_bases,
        contig_cov_sum_bases,
        contig_cov_sum_bases * 100.0 / ref_sum_bases if contig_cov_sum_bases > 0 else 0)
    

def get_gc_and_len_dict(fastafile):
    """Creates a dictionary with the fasta id as key and GC and length as keys
    for the inner dictionary."""
    out_dict = {}
   
    for rec in SeqIO.parse(fastafile, "fasta"):
        out_dict[rec.id] = {}
        out_dict[rec.id]["length"] = len(rec.seq)

    return out_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("coords", help="Output from nucmer show-coords")
    parser.add_argument("reffa", help="Reference fasta file\n")
    parser.add_argument("--cutoff", type=int, default=100, help='Minimum contig length (100)')
    args = parser.parse_args()

    print_metagenome_coverage(args.coords, args.reffa, args.cutoff)
