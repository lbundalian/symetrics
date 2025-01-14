import numpy as np
import pandas as pd
import sqlite3
from sqlite3 import Error
import abc
import logging
from enum import Enum, auto
import json
from sklearn.preprocessing import StandardScaler
from .src.datastruct import *
from .dbcontext import DbContext
from .modelcontext import ModelContext


class ISymetrics(abc.ABC):

    @abc.abstractclassmethod
    def get_silva_score():
        pass

    @abc.abstractclassmethod
    def get_surf_score():
        pass

    @abc.abstractclassmethod
    def get_synvep_score():
        pass

    @abc.abstractclassmethod
    def get_spliceai_score():
        pass

    @abc.abstractclassmethod
    def get_prop_score():
        pass

    @abc.abstractclassmethod
    def get_gnomad_data():
        pass

    @abc.abstractclassmethod
    def get_gnomad_constraints():
        pass

    @abc.abstractclassmethod
    def liftover():
        pass


class Symetrics(ISymetrics):

    _db = None
    _conn = None
    _collection = None
    _gnomad_db = None
    _constraints = None
    _features = ['Synvep', 'SPLICEAI', 'SURF', 'MES', 'GERP', 'CpG', 'CpG_exon', 'RSCU', 'dRSCU', 'F_MRNA', 'F_PREMRNA', 'AF']

    def __init__(self, cfg) -> None:
        with open(cfg, 'r') as file:
            config = json.load(file)

        self._db = DbContext(config['collection']['symetrics'])
        self._gnomad_db = DbContext(config['collection']['gnomad'])
        self._constraints = DbContext(config['collection']['constraints'])
        self._model = ModelContext(config['collection']['model'], self._features)
        self._collection = config

    def execute_query(self, query: str, params: tuple, default_result: dict, columns: list, db: DbContext):
        """
        Executes an SQL query and returns the result as a dictionary.
        If no rows are returned, returns default_result.
        """
        try:
            with db as dbhandler:
                cursor = dbhandler._conn.cursor()
                cursor.execute(query, params)
                row = cursor.fetchone()
                if row:
                    return dict(zip(columns, row))
        except Error as e:
            logging.error(f"Database error: {e}")
        return default_result

    def get_silva_score(self, variant: VariantObject):
        default_values = {
            "CHR": variant._chr,
            "POS": variant._pos,
            "REF": variant._ref,
            "ALT": variant._alt,
            "GENE": "N/A",
            "RSCU": 0.0,
            "dRSCU": 0.0,
            "GERP": 0.0,
            "MES": 0.0,
            "CPG": 0.0,
            "CPGX": 0.0,
            "F_PREMRNA": 0.0,
            "F_MRNA": 0.0
        }
        query = """
            SELECT CHROM, POS, REF, ALT, GENE, GERP, RSCU, dRSCU, CpG, CpG_exon, MES, F_PREMRNA, F_MRNA
            FROM SILVA_SCORE
            WHERE CHROM = ? AND POS = ? AND REF = ? AND ALT = ?
        """
        columns = ["CHR", "POS", "REF", "ALT", "GENE", "GERP", "RSCU", "dRSCU", "CPG", "CPGX", "MES", "F_PREMRNA", "F_MRNA"]
        return self.execute_query(query, (variant._chr, variant._pos, variant._ref, variant._alt), default_values, columns,self._db)

    def get_surf_score(self, variant: VariantObject):
        default_values = {
            "CHR": variant._chr,
            "POS": variant._pos,
            "REF": variant._ref,
            "ALT": variant._alt,
            "SURF": 0.0
        }
        query = """
            SELECT CHR, POS, REF, ALT, SURF
            FROM SURF
            WHERE CHR = ? AND POS = ? AND REF = ? AND ALT = ?
        """
        columns = ["CHR", "POS", "REF", "ALT", "SURF"]
        return self.execute_query(query, (variant._chr, variant._pos, variant._ref, variant._alt), default_values, columns,self._db)

    def get_synvep_score(self, variant: VariantObject):
        default_values = {
            "CHR": variant._chr,
            "POS": variant._pos,
            "REF": variant._ref,
            "ALT": variant._alt,
            "GENE": "N/A",
            "SYNVEP": 0.0
        }
        query = """
            SELECT chr as CHR, pos_GRCh38 as POS, ref as REF, alt as ALT, HGNC_gene_symbol as GENE, synVep as SYNVEP
            FROM SYNVEP
            WHERE chr = ? AND pos_GRCh38 = ? AND ref = ? AND alt = ?
        """
        columns = ["CHR", "POS", "REF", "ALT", "GENE", "SYNVEP"]
        return self.execute_query(query, (variant._chr, variant._pos, variant._ref, variant._alt), default_values, columns,self._db)

    def get_spliceai_score(self, variant: VariantObject):
        default_values = {
            "CHR": variant._chr,
            "POS": variant._pos,
            "REF": variant._ref,
            "ALT": variant._alt,
            "MAX_DS": 0.0
        }
        try:
            with self._db as dbhandler:
                cursor = dbhandler._conn.cursor()
                query = """
                    SELECT chr as CHR, pos as POS, ref as REF, alt as ALT, INFO
                    FROM SPLICEAI
                    WHERE chr = ? AND pos = ? AND ref = ? AND alt = ?
                """
                params = (variant._chr, variant._pos, variant._ref, variant._alt)
                cursor.execute(query, params)
                rows = cursor.fetchall()
                if rows:
                    df = pd.DataFrame(rows, columns=["CHR", "POS", "REF", "ALT", "INFO"])
                    vcf_header = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split('|')
                    df[vcf_header] = df['INFO'].str.split('|', expand=True)
                    df['MAX_DS'] = df[['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']].astype(float).max(axis=1)
                    return df[['CHR', 'POS', "REF", "ALT", "MAX_DS"]].to_dict(orient='records')[0]
        except Error as e:
            logging.error(f"Database error: {e}")
        return default_values

    def get_gnomad_data(self, variant: VariantObject):
        default_values = {
            "CHR": variant._chr,
            "POS": variant._pos,
            "REF": variant._ref,
            "ALT": variant._alt,
            "AC": 0,
            "AN": 0,
            "AF": 0.0
        }
        query = """
            SELECT chr as CHR, pos as POS, ref as REF, alt as ALT, AC, AN, AF
            FROM gnomad_db
            WHERE chr = ? AND pos = ? AND ref = ? AND alt = ?
        """
        columns = ["CHR", "POS", "REF", "ALT", "AC", "AN", "AF"]
        return self.execute_query(query, (variant._chr, variant._pos, variant._ref, variant._alt), default_values, columns, self._gnomad_db)

    def get_gnomad_constraints(self, gene=''):
        default_values = {
            "gene": gene,
            "transcript": "N/A",
            "syn_z": 0.0,
            "mis_z": 0.0,
            "lof_z": 0.0,
            "pLI": 0.0
        }
        query = """
            SELECT * FROM GNOMADv4Constraints WHERE gene = ?
        """
        columns = ["gene", "transcript", "syn_z", "mis_z", "lof_z", "pLI"]
        return self.execute_query(query, (gene,), default_values, columns, self._constraints)

    def liftover(self, variant: VariantObject):
        """
        Perform liftover depending on the assembly.
        If the assembly is hg19, return hg38; otherwise, return hg19.
        """
        liftover_variant = None
        try:
            with self._db as dbhandler:
                cursor = dbhandler._conn.cursor()
                if variant._genome == GenomeReference.hg19:
                    query = """
                        SELECT chr as CHR, pos_GRCh38 as POS, ref as REF, alt as ALT
                        FROM SYNVEP
                        WHERE chr = ? AND pos = ? AND ref = ? AND alt = ?
                    """
                    new_reference = GenomeReference.hg38
                else:
                    query = """
                        SELECT chr as CHR, pos as POS, ref as REF, alt as ALT
                        FROM SYNVEP
                        WHERE chr = ? AND pos_GRCh38 = ? AND ref = ? AND alt = ?
                    """
                    new_reference = GenomeReference.hg19
                params = (variant._chr, variant._pos, variant._ref, variant._alt)
                cursor.execute(query, params)
                row = cursor.fetchone()
                if row:
                    liftover_variant = VariantObject(
                        chr=row[0],
                        pos=row[1],
                        ref=row[2],
                        alt=row[3],
                        genome=new_reference
                    )
        except Error as e:
            logging.error(f"Liftover error: {e}")
        return liftover_variant or variant


    def get_prop_score(self, group=MetricsGroup.SYNVEP.name, gene=''):
        default_values = {
            "GENE": gene,
            "PVAL": 1.0,
            "FDR": 1.0,
            "SYMETRIC_SCORE": 0.0,
            "NORM_SYMETRIC_SCORE": 0.0
        }
        if group in MetricsGroup.__members__:
            query = f"""
                SELECT GENES as GENE, pval as PVAL, fdr as FDR, z as SYMETRIC_SCORE, norm_z as NORM_SYMETRIC_SCORE
                FROM GNOMADv4{group}
                WHERE GENES = ?
            """
            columns = ["GENE", "PVAL", "FDR", "SYMETRIC_SCORE", "NORM_SYMETRIC_SCORE"]
            return self.execute_query(query, (gene,), default_values, columns,self._constraints)
        else:
            logging.error(f"Group: {group} is skipped")
        return default_values

    def get_all_prop_score(self, group=MetricsGroup.SYNVEP.name):
        scores = (0, 1)
        if group in MetricsGroup.__members__:
            try:
                with self._constraints as dbhandler:
                    cursor = dbhandler._conn.cursor()
                    query = f"""
                        SELECT z as SYMETRIC_SCORE
                        FROM GNOMADv4{group}
                    """
                    cursor.execute(query)
                    rows = cursor.fetchall()
                    rows = [r[0] for r in rows]
                    scores = (np.mean(rows), np.std(rows))
            except Error as e:
                logging.error(f"Database error: {e}")
        else:
            logging.error(f"Group: {group} is skipped")
        return scores

    def get_variant_list(self, gene: str):
        """
        Get a list of variants associated with the gene.
        """
        variant_list = []
        try:
            with self._db as dbhandler:
                cursor = dbhandler._conn.cursor()
                query = """
                    SELECT chr as CHR, pos_GRCh38 as POS, pos as POS_HG19, ref as REF, alt as ALT, HGNC_gene_symbol as GENE, synVep as SYNVEP
                    FROM SYNVEP
                    WHERE HGNC_gene_symbol = ?
                """
                params = (gene,)
                cursor.execute(query, params)
                rows = cursor.fetchall()
                if rows:
                    columns = ["CHR", "POS", "POS_HG19", "REF", "ALT", "GENE", "SYNVEP"]
                    variant_list = [dict(zip(columns, row)) for row in rows]
        except Error as e:
            logging.error(f"Database error: {e}")
        return variant_list

    def predict_probability(self, scores: dict):
        pred = None
        try:
            pred = self._model.predict(scores)
        except Exception as e:
            logging.error(f"Model error: {e}")
        return pred

    def get_variant_batch(self, gene: str, start: int, end: int):
        """
        Get a batch of variants associated with the gene.
        """
        variant_list = []
        try:
            with self._db as dbhandler:
                cursor = dbhandler._conn.cursor()
                query = """
                    SELECT chr as CHR, pos_GRCh38 as POS, pos as POS_HG19, ref as REF, alt as ALT, HGNC_gene_symbol as GENE, synVep as SYNVEP
                    FROM SYNVEP
                    WHERE HGNC_gene_symbol = ?
                    LIMIT ? OFFSET ?
                """
                params = (gene, end - start, start)
                cursor.execute(query, params)
                rows = cursor.fetchall()
                if rows:
                    columns = ["CHR", "POS", "POS_HG19", "REF", "ALT", "GENE", "SYNVEP"]
                    variant_list = [dict(zip(columns, row)) for row in rows]
        except Error as e:
            logging.error(f"Database error: {e}")
        return variant_list
