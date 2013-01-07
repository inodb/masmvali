import numpy as np


class EqualNameImpliesEquality:
    """Inheriting from EqualNameImpliesEquality allows using instances of the
    class in a set. If the name attribute of the instance is equal to the name
    of the other instance they are considered equal."""
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.name == other.name
        else:
            raise(Exception('Comparing different classes on name.'))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.name)


class Reference(EqualNameImpliesEquality):
    """Reference keeps track of the phyologeny of a reference based on NCBI's taxonomy."""
    TAX_LVLS = ['strain', 'species', 'genus', 'family', 'order', 'class',
                'phylum', 'superphylum', 'superkingdom', 'life']

    def __init__(self, phylogeny):
        """Expects an ordered list of phylogeny from lowest rank to highest as TAX_LVLS."""
        assert(len(phylogeny) == len(self.TAX_LVLS))

        # Could have been a named tuple, but 'class' is not allowed as attribute name
        self.phyl = dict(zip(self.TAX_LVLS, phylogeny))
        self.name = self.phyl["strain"]

    def get_highest_lca(self, other_strains, missing_value="N/A"):
        """Get the highest taxonomic rank of the Lowest Common Ancestor between
        other_strains and self."""
        tax_order = dict(zip(self.TAX_LVLS, range(len(self.TAX_LVLS))))
        max_tl = 0
        for os in other_strains:
            max_tl = max(tax_order[self.get_lca(os, missing_value=missing_value)[0]], max_tl)

        return((self.TAX_LVLS[max_tl], self.phyl[self.TAX_LVLS[max_tl]]))

    def get_lca(self, strain, missing_value="N/A"):
        """Get the Lowest Common Ancestor between strain and self.

        strain -- a Reference instance
        """
        assert(isinstance(strain, self.__class__))

        for rank in self.TAX_LVLS:
            if self.phyl[rank] == strain.phyl[rank] != missing_value:
                return (rank, strain.phyl[rank])

        raise(Exception("No LCA? Not a cellular organism?"))

    def __str__(self):
        return(self.name)


class ReferenceSet():
    """A set of references"""
    def __init__(self, refphylfile, refstatsfile, missing_value="N/A"):
        # Get reference length/GC info
        reflens = readtable(refstatsfile, sep="\t")

        # Get phylogeny info
        rp = np.genfromtxt(refphylfile, names=True, dtype=None,
                           missing_values=missing_value, delimiter="\t")
        tax_lvls = [rp.dtype.names[i] for i in range(4, 14)]

        # Create reference dic
        self.refs = {}
        for refrec in rp:
            # Only use references that actually have fasta sequences
            fn = refrec["fasta_name"]
            if fn != missing_value:
                # Phylogeny info
                phylogeny = [refrec[t] for t in reversed(tax_lvls)]
                r = Reference(phylogeny)

                # Reference length/GC/ratio covered
                r.length = reflens[fn]["length"]
                r.ratio_covered = reflens[fn]["ratio_covered"]
                r.cov = reflens[fn]["cov"]
                r.gc_content = reflens[fn]["GC_content"]

                # refs indexable by genome id, fasta name and reference name
                self.refs[refrec["fasta_name"]] = r
                self.refs[refrec["gen_id"]] = r
                self.refs[r.name] = r

    def get(self, name):
        return self.refs[name]

    def __len__(self):
        return len(set(self.refs.values()))


def readtable(tablefile, sep=None):
    """Reads given table separated 'sep'.

    Returns a two-dimensional dictionary. Outer dictionary uses first column in
    'tablefile' as key, inner dictionary uses column names found in the first
    line as column names. The number of column names should be equal to number
    of columns - 1."""
    table2d = {}

    # Get column names
    tfh = open(tablefile, "r")
    line = tfh.readline()
    cols = line.strip().split(sep)

    # Insert rows
    for line in tfh:
        splits = line.strip().split(sep)
        table2d[splits[0]] = {}
        for i in range(1, len(splits)):
            table2d[splits[0]][cols[i - 1]] = splits[i]

    return(table2d)
