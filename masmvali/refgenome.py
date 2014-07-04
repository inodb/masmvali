import numpy as np
from utils import readtable


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
    def __init__(self, refphylfile, refstatsfile, contigs_to_refs_dict, missing_value="N/A"):
        # Get reference length/GC info
        reflens = readtable(refstatsfile, sep="\t")

        # Get phylogeny info
        rp = np.genfromtxt(refphylfile, names=True, dtype=None,
                           missing_values=missing_value, delimiter="\t")
        tax_lvls = [rp.dtype.names[i] for i in range(1, 11)]

        # Create reference dic
        self.refs = {}
        refs_in_contigs_to_refs_dict = set(contigs_to_refs_dict.itervalues())
        for refrec in rp:
            # Only use references that actually have fasta sequences
            if refrec["topname"] in refs_in_contigs_to_refs_dict:
                # Phylogeny info
                phylogeny = [refrec[t] for t in reversed(tax_lvls)]
                r = Reference(phylogeny)

                r.length = 0
                r.ratio_covered = 0
                r.cov = 0
                r.gc_content = 0
                for c in contigs_to_refs_dict:
                    if contigs_to_refs_dict[c] == refrec["topname"]:
                        contig_length = int(reflens[c]["length"])
                        r.length += contig_length
                        #print c, reflens[c]
                        #TODO: the refstats stuff should be changed to
                        # actually take a reference fasta and a bam file
                        r.ratio_covered += float(reflens[c]["ratio_covered"]) * contig_length
                        r.cov += float(reflens[c]["cov"]) * contig_length
                        r.gc_content += float(reflens[c]["GC_content"]) * contig_length
                if r.length == 0:
                    raise(Exception("No contigs found for %s in contigs_to_refs_dict" % refrec["topname"]))
                # Divide by sum of contig lengths belonging to reference
                r.ratio_covered = r.ratio_covered / r.length
                r.cov = r.cov / r.length
                r.gc_content = r.cov / r.length

                # refs indexed by reference name
                self.refs[r.name] = r

        self.mg_length = sum([r.length for r in self])

    def get(self, name):
        return self.refs[name]

    def __len__(self):
        return len(set(self.refs.itervalues()))

    def __iter__(self):
        return iter(set(self.refs.itervalues()))
