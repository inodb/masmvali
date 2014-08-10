import pandas as pd
import numpy as np
from common.cache import property_cached

def strip(text):
    try:
        return text.strip()
    except AttributeError:
        return text


def list_it(text):
    return [(int(c[0]), int(c[1])) for c in [s.split(':') for s in text.split()]]


class Taxonomy(object):
    def __init__(self, node_file, taxonomy_file, ref_taxonomy_ids_file):
        self.nodes = pd.read_csv(node_file, sep="|",
                                 skipinitialspace=True, index_col=False,
                                 names=["tax_id","parent_tax_id","rank","embl_code",
                                        "division_id","inherited_div_flag","genetic_code_id",
                                        "inherited_GC_flag", "mitochondrial genetic_code_id",
                                        "inherited_MGC_flag", "Genbank_hidden_flag","hidden_subtree_flag",
                                        "comments"],
                                converters={"rank":strip,"embl_code":strip}).set_index("tax_id")
        self.tax_names = pd.read_csv(taxonomy_file, sep="|",
                                     skipinitialspace=True, index_col=False,
                                     names=["tax_id","name_txt","uniq_name","name_class"],
                                     converters={"name_txt":strip,"uniq_name":strip,"name_class":strip})
        self.tids = [int(t.rstrip('\n')) for t in open(ref_taxonomy_ids_file).readlines()]

    def get_rank(self, tax_id):
        """Get the rank of a given tax-id. Get the parent rank recursively if it has no rank."""
        if tax_id == 1:
            rank = "domain"
        elif tax_id == 0:
            rank = "non-existent"
        else:
            node_entry = self.nodes.ix[tax_id]
            rank = node_entry["rank"]
            if rank == "no rank":
                parent_rank = self.nodes.ix[node_entry.parent_tax_id]["rank"]
                if parent_rank == "no rank":
                    rank = self.get_rank(node_entry["parent_tax_id"])
                else:
                    rank = parent_rank
        return rank

    def get_full_lineage(self, tid):
        """Traverses all parent taxonomy ids to return a list from lowest taxonomy id to highest"""
        ptid = self.nodes.ix[tid]["parent_tax_id"]
        if ptid == tid:
            return [tid]
        else:
            return [tid] + self.get_full_lineage(ptid)

    def get_lca(self, lin, linref):
        """Searches taxonomy id LCA of lin in linref. Lin and linref are expected to be ordered from lowest taxonomy id to highest."""
        for tid in lin:
            try:
                i = linref.index(tid)
                return linref[i]
            except ValueError:
                continue
        raise(Exception("No LCA found? Should always be Life!"))

    def get_lca_rank(self, row):
        """Get the rank of each LCA for all kmers that are not the same as the one the contig has been classified as."""
        linref = self.get_full_lineage(row.name)
        rv = pd.Series(["unknown"] * len(row), index=row.index)
        for i in row.index:
            lin = self.get_full_lineage(i)
            lca_taxid = self.get_lca(lin, linref)
            if i == lca_taxid:
                # prefix the rank of a taxid that is contained in the lineage of linref as 'contained_'
                rank = "contained_" + self.get_rank(lca_taxid)
            else:
                rank = self.get_rank(lca_taxid)
            rv[i] = rank
        return rv

    def lca_rank_count_row(self, row, unique_ranks, taxa2lca_rank):
        """Count number of kmers for each LCA rank from kraken dataframe"""
        rv = pd.Series([0]*(len(unique_ranks)+2),index=list(unique_ranks) + ["correct","unknown"])
        for (tid, count) in row.lca:
            if not row.classified:
                raise(Exception("Unclassified contigs are not implemented"))
            if tid == 0:
                rv["unknown"] += count
            elif row.tax_id != tid:
                rank = taxa2lca_rank.ix[row.tax_id, [tid]].values[0]
                rv[rank] += count
            else:
                rv["correct"] += count

        return rv

    def get_lca_rank_kmer_count_sort(self, krakenpd):
        get_uniq_tids = lambda x: pd.Series(index=set([tid for (tid, count) in x]))
        uniq_tids = list(krakenpd[krakenpd.length > 500].lca.apply(get_uniq_tids).columns)
        taxa_lca_rank = pd.DataFrame(index=sorted(self.tids),
                                     columns=sorted([t for t in uniq_tids if t != 0])
                                    ).fillna(value="unknown")
        taxa2lca_rank = taxa_lca_rank.apply(self.get_lca_rank, axis=1)
        unique_ranks = np.unique(taxa2lca_rank.as_matrix())
        lca_rank_kmer_count = krakenpd[(krakenpd.length > 500) & krakenpd.tax_id.isin(self.tids)].apply(self.lca_rank_count_row,args=(unique_ranks, taxa2lca_rank),axis=1).fillna(0)
        lca_rank_kmer_count_sort = lca_rank_kmer_count.reindex_axis(reversed(['unknown'] +
              ["contained_" + r for r in ['domain','superkingdom','phylum','class','order','family','genus','species']] +
              ['domain','superkingdom','superphylum','phylum','class','order','family','genus','species']), axis=1)

        return lca_rank_kmer_count_sort


class Kraken(object):
    def __init__(self, kraken_file, taxonomy):
        # Taxonomy object
        self.taxonomy = taxonomy
        self.df = pd.read_csv(kraken_file, sep="\t",
                              true_values=["C"], false_values=["U"],
                              names=["classified","seq_id","tax_id","length","lca"],
                              converters={"lca":list_it}).sort("seq_id").reset_index().set_index("seq_id")

    @property_cached
    def lca_rank_kmer_count_sort(self):
        return self.taxonomy.get_lca_rank_kmer_count_sort(self.df)
