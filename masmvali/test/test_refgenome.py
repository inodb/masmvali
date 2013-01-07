from refgenome import Reference, ReferenceSet
from utils import get_shell_output, get_testfile, make_dir, get_outdir


def test_referenceset():
    refstatsfile = get_testfile('cm-500pgun-ref-stats')
    refphylfile = get_testfile('cm-ref-phyl')

    r = Reference(['Acidobacterium capsulatum ATCC 51196', 'Acidobacterium capsulatum', 'Acidobacterium', 'Acidobacteriaceae', 'Acidobacteriales', 'Acidobacteriia', 'Acidobacteria', 'Fibrobacteres/Acidobacteria group', 'Bacteria', 'cellular organisms'])
    assert(r.phyl["life"] == "cellular organisms")

    refs = ReferenceSet(refphylfile, refstatsfile)
    # There are 59 references, all indexable by id, fasta name or strain name
    assert(len(refs) == 59)
    assert(len(refs.refs) == 177)

    # Every reference's full phylogeny is in the phyl dictionary attribute
    assert(refs.get("Acidobacterium capsulatum ATCC 51196").phyl["life"] == "cellular organisms")

    # The LCA of Shewanella baltica strains should be on the species level
    sh1 = refs.get("Shewanella_baltica_OS185")
    sh2 = refs.get("Shewanella_baltica_OS223,")
    assert(sh1.get_lca(sh2) == ('species', 'Shewanella baltica'))

    # LCA of an Archaea and a Bacteria should be life
    ar = refs.get("Archaeoglobus_fulgidus_DSM_4304")
    th = refs.get("gi|222528057|ref|NC_012034.1|")
    assert(ar.get_lca(th) == ('life', 'cellular organisms'))

    # Highest LCA of Bacteria and Archaea is Life
    assert(ar.get_highest_lca([th, sh1, sh2]) == ('life', 'cellular organisms'))
    # Highest LCA of given bacteria is domain or superkingdom in NCBI naming
    assert(th.get_highest_lca([sh1, sh2]) == ('superkingdom', 'Bacteria'))


def main():
    test_referenceset()

if __name__ == "__main__":
    main()
