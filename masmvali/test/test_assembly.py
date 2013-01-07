from assembly import get_read_contig_mappings, AssemblyValidation
from refgenome import ReferenceSet
from utils import get_shell_output, get_testfile, make_dir, get_outdir


def test_read_contig_mappings():
    bamasm = get_testfile('cm-500pgun-asm-b2mv31-bam')
    contigfa = get_testfile('cm-500pgun-asm-b2mv31-fa')
    bamref = get_testfile('cm-500pgun-ref-bam')
    refstatsfile = get_testfile('cm-500pgun-ref-stats')
    refphylfile = get_testfile('cm-ref-phyl')

    refs = ReferenceSet(refphylfile, refstatsfile)
    reads, contigs = get_read_contig_mappings(bamref, bamasm, refs, contigfa)


def test_assemblyvalidation():
    bamasm = get_testfile('cm-500pgun-asm-b2mv31-bam')
    contigfa = get_testfile('cm-500pgun-asm-b2mv31-fa')
    bamref = get_testfile('cm-500pgun-ref-bam')
    refstatsfile = get_testfile('cm-500pgun-ref-stats')
    refphylfile = get_testfile('cm-ref-phyl')
    nucmercoords = get_testfile('cm-500pgun-val-nucmer')

    val = AssemblyValidation(bamref, bamasm, refphylfile, refstatsfile, contigfa, nucmercoords)
    assert(len(val.contigs) == int(get_shell_output("grep -c '^>' " + contigfa)[0]))

    make_dir(get_outdir() + "masm")
    with open(get_outdir() + "masm" + "/contig-purity.tsv", "w") as fh:
        fh.write(val.contigs.itervalues().next().str_stats_header() + "\n")
        for c in val.contigs.itervalues():
            fh.write(c.str_stats() + "\n")


def main():
    #test_read_contig_mappings()
    test_assemblyvalidation()


if __name__ == "__main__":
    main()
