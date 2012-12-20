from refgenome import Reference, ReferenceSet


def test_referenceset():
    refphylfile = '/bubo/glob/g16/inod/metagenomics/result/chris-mock/phylogeny-references.tsv'
    refstatsfile = '/bubo/glob/g16/inod/metagenomics/result/chris-mock/Project_ID793_dAmore/Sample_500pg_unbalanced/ma2-out/reference-stats/ref.stats'

    r = Reference(['Acidobacterium capsulatum ATCC 51196', 'Acidobacterium capsulatum', 'Acidobacterium', 'Acidobacteriaceae', 'Acidobacteriales', 'Acidobacteriia', 'Acidobacteria', 'Fibrobacteres/Acidobacteria group', 'Bacteria', 'cellular organisms'])
    assert(r.phyl["life"] == "cellular organisms")

    refs = ReferenceSet(refphylfile, refstatsfile)
    assert(len(refs) == 59)
    assert(len(refs.refs) == 177)
    assert(refs.get("Acidobacterium capsulatum ATCC 51196").phyl["life"] == "cellular organisms")


def main():
    test_referenceset()

if __name__ == "__main__":
    main()
