import sqlite3

class Coords:
    def __init__(self, coordsfile):
        self.coordsfile = coordsfile
        self.cursor = get_coords_db_cursor(coordsfile)

    def calc_genome_contig_cov_in_bases(self, cut_off=100):
        return(calc_genome_contig_cov_in_bases(self.cursor, cut_off))

    def calc_max_aln_purity_per_contig(self, contigs=None, cut_off=100):
        return(calc_max_aln_purity_per_contig(self.cursor, contigs=contigs, cut_off=cut_off))

    def calc_alignedbases_per_contig(self, contigs=None, cut_off=100):
        return(calc_alignedbases_per_contig(self.cursor, contigs=contigs, cut_off=cut_off))


def get_coords_db_cursor(coordsfile):
    # Create db
    #dbc = sqlite3.connect(':memory:')
    dbc = sqlite3.connect(coordsfile + ".sqlite")
    try:
        dbc.execute("""Create table Coords (ID INTEGER PRIMARY KEY, S1 INTEGER, E1
                    INTEGER, S2 INTEGER, E2 INTEGER, LEN1 INTEGER, LEN2 INTEGER,
                    IDY REAL, LENR INTEGER, LENQ INTEGER, COVR REAL, COVQ REAL,
                    REFID TEXT, QRYID TEXT)""")
        dbc.row_factory = sqlite3.Row

        # Parse file and add to db
        columns = ('S1', 'E1', 'S2', 'E2', 'LEN1', 'LEN2', 'IDY', 'LENR',
                   'LENQ', 'COVR', 'COVQ', 'REFID', 'QRYID')
        cfh = open(coordsfile, "r")
        for line in cfh:
            parsed_line = dict(zip(columns, line.split()))
            dbc.execute("""Insert into Coords values(NULL, :S1, :E1, :S2, :E2,
                        :LEN1, :LEN2, :IDY, :LENR, :LENQ, :COVR, :COVQ, :REFID,
                        :QRYID)""", parsed_line)
    except:
        dbc.row_factory = sqlite3.Row
        pass
    dbc.commit()

    return dbc.cursor()

def calc_genome_contig_cov_in_bases(cursor, cut_off=100):
    """Genome contig coverage is a metric that indicates how well the genome is
    covered by contigs. For each contig only the purest alignment is
    considered. Purity is defined as COVQ * IDY / 10,000. If there are muliple
    alignments for a contig with maximum purity, both are used. Covered bases
    are only counted once and computed by multiplying the length of the
    alignment in the reference (S1 - E1 + 1) by the alignment identity (IDY).
    If contigs overlap the first contig based on its location in the reference
    genome is used. If the second contig extends further than the first, the
    second's bases are added as well using the IDY of the second. One might
    prefer to use the purest alignment in case of overlap but I'll implement
    that only if persuaded with dinner.

    Returns a dictionary of the genome contig coverage per genome as the number
    of non-overlapping bases that are covered by contigs aligning with maximum
    purity.
    """
    refcov = {}  # nr of bases covered by contigs aligned with max purity
    prev_e1 = 0
    for row in cursor.execute("""SELECT * FROM (SELECT *, max(COVQ * IDY /
                              10000) AS purity FROM Coords GROUP BY QRYID)
                              WHERE COVQ * IDY / 10000 == purity AND LENQ >= :cut_off ORDER BY
                              REFID, S1 ASC""", dict(cut_off=cut_off)):
        if row["REFID"] in refcov:
            if prev_e1 >= row["E1"]:
                continue
            elif prev_e1 >= row["S1"]:
                refcov[row["REFID"]] += int((row["E1"] -
                                            prev_e1) * (row["IDY"] / 100.0))
            else:
                refcov[row["REFID"]] += int((row["E1"] - row["S1"] + 1) *
                                            (row["IDY"] / 100.0))
        else:
            refcov[row["REFID"]] = int((row["E1"] - row["S1"])
                                       * (row["IDY"] / 100.0)) + 1
        prev_e1 = row["E1"]

    return(refcov)


def calc_max_aln_purity_per_contig(cursor, contigs=None, cut_off=100):
    query = """SELECT * FROM (SELECT *, max(COVQ * IDY / 10000) AS purity FROM
            Coords GROUP BY QRYID) WHERE COVQ * IDY / 10000 == purity AND LENQ
            >= :cut_off"""

    if contigs is None:
        q_max_aln_purity = {}
        for row in cursor.execute(query, dict(cut_off=cut_off)):
            q_max_aln_purity[row["QRYID"]] = dict(max_aln_purity=row["purity"],
                                                length=row["LENQ"])
    else:
        q_max_aln_purity = contigs
        for row in cursor.execute(query, dict(cut_off=cut_off)):
            # Length of the contig should already be set if this is a ContigSet
            # by parsing the fasta file
            # TODO: assertion fails because nucmer includes Ns
            #assert(contigs[row["QRYID"]].length == row["LENQ"])

            q_max_aln_purity[row["QRYID"]].max_aln_purity = row["purity"]
            q_max_aln_purity[row["QRYID"]].max_aln_strain = row["REFID"]

    return(q_max_aln_purity)


def calc_alignedbases_per_contig(cursor, contigs=None, cut_off=100):
    query = """SELECT *, min(S2, E2) AS start, max(S2, E2) as end FROM COORDS
            WHERE LENQ >= :cut_off ORDER BY QRYID, start ASC"""

    if contigs is None:
        q_aln_bases = {}  # nr of bases aligned per contig
        prev_end = 0
        for row in cursor.execute(query, dict(cut_off=cut_off)):

            if not row["QRYID"] in q_aln_bases:
                q_aln_bases[row["QRYID"]] = row["end"] - row["start"]
            else:
                if prev_end >= row["end"]:
                    continue
                elif prev_end >= row["start"]:
                    q_aln_bases[row["QRYID"]] += row["end"] - prev_end
                else:
                    q_aln_bases[row["QRYID"]] += row["end"] - row["start"] + 1
            prev_end = row["end"]

        return(q_aln_bases)
    else:
        assert(sum([hasattr(c, "aln_bases") for c in contigs]) == 0)

        prev_end = 0
        for row in cursor.execute(query, dict(cut_off=cut_off)):
            if not hasattr(contigs[row["QRYID"]], "aln_bases"):
                contigs[row["QRYID"]].aln_bases = row["end"] - row["start"]
            else:
                if prev_end >= row["end"]:
                    continue
                elif prev_end >= row["start"]:
                    contigs[row["QRYID"]].aln_bases += row["end"] - prev_end
                else:
                    contigs[row["QRYID"]].aln_bases += row["end"] - row["start"] + 1
            prev_end = row["end"]

        return(contigs)
