import numpy as np
import pandas as pd
import sqlite3
from sqlite3 import Error
import abc
import logging
from enum import Enum,auto
import json
from sklearn.preprocessing import StandardScaler
from .src.datastruct import *
from .dbcontext import DbContext


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
    _collection = None
    
    def __init__(self, cfg) -> None:

        with open(cfg, 'r') as file:
            config = json.load(file)
    
        self._db = DbContext(config['collection']['symetrics'])
        self._gnomad_db = DbContext(config['collection']['gnomad'])
        self._constraints = DbContext(config['collection']['constraints'])
        self._collection = config

    
    def get_silva_score(self,variant: VariantObject):

        """
        
        Get the RSCU, dRSCU, GERP and CpG/CpG_Exon of a given variant (reference: hg19)
        
        Args:
            variant: A VariantObject instance representing the chromosome, position, reference allele and alternative allele of a variant
        
        Returns:
            silva_scores: A dictionary returning the scores along with the variant information.

        Examples:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> variant = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg19)
            >>> silva = symetrics.get_silva_score(variant)
        
        """

        silva_scores = None
        try:
            # dont forget silva is hg19
            with self._db as dbhandler:
                
                silva_cursor = dbhandler._conn.cursor()
                
                silva_query = """
                        SELECT * FROM SILVA_SCORE 
                        WHERE CHROM = ? AND POS = ? AND REF = ? AND ALT = ?
                        """
                params = (variant._chr, variant._pos, variant._ref, variant._alt)

                silva_cursor.execute(silva_query,params)
                silva_rows = silva_cursor.fetchall()
                silva_scores = silva_rows[0]
                silva_scores = {
                    "CHR": silva_scores[0],
                    "POS": silva_scores[1],
                    "REF": silva_scores[3],
                    "ALT": silva_scores[4],
                    "GENE": silva_scores[5],
                    "RSCU": silva_scores[8],
                    "dRSCU": silva_scores[9],
                    "GERP": silva_scores[7],
                    "MES": silva_scores[12],
                    "CPG": silva_scores[10],
                    "CPGX": silva_scores[11]
                }

        except Error as e:
            logging.error(f"Connection to {self._db} failed")
        


        return silva_scores

    def get_surf_score(self,variant: VariantObject):
        
        """
        
        Get the SURF a given variant (reference: hg38)
        
        Args:
            variant: A VariantObject instance representing the chromosome, position, reference allele and alternative allele of a variant
        
        Returns:
            surf_scores: A dictionary returning the scores along with the variant information.
        
        Examples:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> variant = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg38)
            >>> surf = symetrics.get_surf_score(variant)
        

        """

        surf_scores = None
        try:
            # SURF is hg38
            with self._db as dbhandler:
                
                surf_cursor = dbhandler._conn.cursor()
                
                surf_query = """
                    SELECT CHR, POS, REF, ALT, GENE, SURF
                    FROM SURF
                    WHERE CHR = ? AND POS = ? AND REF = ? AND ALT = ?
                    """
                
                params = (variant._chr, variant._pos, variant._ref, variant._alt)

                surf_cursor.execute(surf_query,params)
                surf_rows = surf_cursor.fetchall()
                surf_scores = surf_rows[0]
                surf_scores = {
                    "CHR": surf_scores[0],
                    "POS": surf_scores[1],
                    "REF": surf_scores[2],
                    "ALT": surf_scores[3],
                    "SURF": surf_scores[5]

                }
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return surf_scores
    
    def get_synvep_score(self,variant: VariantObject):

        """
        
        Get the SYNVEP a given variant (reference: hg38/hg19)
        https://services.bromberglab.org/synvep/home
        
        Args:
            variant: A VariantObject instance representing the chromosome, position, reference allele and alternative allele of a variant
        
        Returns:
            synvep_scores: A dictionary returning the scores along with the variant information.
        
        Examples:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> variant_hg19 = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg19)
            >>> variant_hg38 = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg38)
            >>> synvep_hg19 = symetrics.get_synvep_score(variant_hg19)
            >>> synvep_hg38 = symetrics.get_synvep_score(variant_hg38)

        """

        synvep_scores = None

        try:
            # synvep is hg38 (pos_GRCh38) abd hg19 (pos)
            with self._db as dbhandler:
                synvep_cursor = dbhandler._conn.cursor()
                synvep_query = ''
                if variant._genome.name == GenomeReference.hg38.name:
                    synvep_query = """
                        SELECT chr as CHR, pos_GRCh38 as POS, ref as REF, alt as ALT, HGNC_gene_symbol as GENE, synVep as SYNVEP
                        FROM SYNVEP
                        WHERE chr = ? AND pos_GRCh38 = ? AND ref = ? AND alt = ?
                        """    
                    
                elif variant._genome.name == GenomeReference.hg19.name:                    
                    synvep_query = """
                        SELECT chr as CHR, pos_GRCh38 as POS, ref as REF, alt as ALT, HGNC_gene_symbol as GENE, synVep as SYNVEP
                        FROM SYNVEP
                        WHERE chr = ? AND pos = ? AND ref = ? AND alt = ?
                        """
                params = (variant._chr, variant._pos, variant._ref, variant._alt)
                
                synvep_cursor.execute(synvep_query,params)
                synvep_rows = synvep_cursor.fetchall()
                synvep_scores = synvep_rows[0]
                synvep_scores = {
                    "CHR": synvep_scores[0],
                    "POS": synvep_scores[1],
                    "REF": synvep_scores[2],
                    "ALT": synvep_scores[3],
                    "GENE": synvep_scores[4],
                    "SYNVEP": synvep_scores[5]

                }
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return synvep_scores

    def get_spliceai_score(self, variant: VariantObject):

        """
        
        Get the SpliceAI a given variant (reference: hg38)
        https://spliceailookup.broadinstitute.org/
        
        Args:
            variant: A VariantObject instance representing the chromosome, position, reference allele and alternative allele of a variant
        
        Returns:
            spliceai_score: A dictionary returning the scores along with the variant information.
        
        Examples:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> variant = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg38)
            >>> spliceai = symetrics.get_spliceai_score(variant)

        """
                        
        spliceai_score = None
        try:
            # synvep is hg38 (pos_GRCh38) abd hg19 (pos)
            with self._db as dbhandler:
                spliceai_cursor = dbhandler._conn.cursor()
                                
                spliceai_query = """
                SELECT chr as CHR, pos as POS, ref as REF, alt as ALT, INFO
                FROM SPLICEAI
                WHERE chr = ? AND pos = ? AND ref = ? AND alt = ?
                """

                # Prepare the parameters tuple
                params = (variant._chr, variant._pos, variant._ref, variant._alt)

                spliceai_cursor.execute(spliceai_query,params)
                spliceai_rows = spliceai_cursor.fetchall()
                spliceai_score = pd.DataFrame(spliceai_rows)
                spliceai_score.columns = ['CHR','POS','REF','ALT','INFO']
                if not spliceai_score.empty:
                    vcf_header = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"
                    vcf_header = vcf_header.split('|')
                    spliceai_score[vcf_header] = spliceai_score['INFO'].str.split('|', expand=True)
                    spliceai_score['MAX_DS'] = spliceai_score.apply(lambda row: max(row['DS_AG'],row['DS_AL'],row['DS_DG'],row['DS_DL']), axis=1)
                    spliceai_score = spliceai_score[['CHR','POS','REF','ALT','MAX_DS']]
                    spliceai_score = spliceai_score.to_dict(orient='records')

        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return spliceai_score

    def get_prop_score(self,group = MetricsGroup.SYNVEP.name,gene = ''):
        
        """
        
        Get the SYMETRICS score for a given gene abd metrics group. The score was calculated from the pooled z proportion test of different
        metrics group with their corresponding threhold:
            - SYNVEP: 0.5
            - GERP: 4
            - CpG: 
            - CpG_exon: 1
            - RSCU:
            - dRSCU:
            - SpliceAI: 0.8
            - SURF: 0.3
        
        Args:
            gene: A string representing the HGNC Symbol of a gene
        
        Returns:
            scores: A dictionary returning the pvalues and fdr acquired from the test and the score before and after scaling.
        
        Examples:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> score = symetrics.get_prop_score(group = 'SYNVEP',gene = 'A1BG')
        
        """

        scores = None
        scaler = StandardScaler()
        

        if group in MetricsGroup.__members__:
            scores = None
            try:
                with self._constraints as dbhandler:
                            
                    cursor = dbhandler._conn.cursor()
                            
                    query = f"SELECT GENES as GENE, pval as PVAL, fdr as FDR, z as SYMETRIC_SCORE, norm_z as NORM_SYMETRIC_SCORE FROM GNOMADv4{group} WHERE GENES = ?"
                            
                    params = (gene,)

                    cursor.execute(query,params)
                    rows = cursor.fetchall()
                    scores = rows[0]
                    scores = {
                                "GENE": scores[0],
                                "PVAL": scores[1],
                                "FDR": scores[2],
                                "SYMERIC_SCORE": scores[3],
                                "NORM_SYMERIC_SCORE": scores[4]
                    }

            except Error as e:
                    logging.error(e)
                    logging.error(f"Connection to {self._constraints} failed")
                
                    
        else:
            logging.error(f'Group: {group} is not valid')       
    
        return scores    
    
    def get_gnomad_data(self, variant: VariantObject):

        """
        
        Get the gnomad information related to the alleles of the given variant (allele count, allele number and allele frequency)
        
        Args:
            variant: A VariantObject instance representing the chromosome, position, reference allele and alternative allele of a variant
        
        Returns:
            gnomad_data: A dictionary containing the AC, AN, AF and variant information
        
        Examples:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> variant_hg38 = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg38)
            >>> gnomad_hg38 = symetrics.get_gnomad_data(variant_hg38)

        """


        gnomad_conn = None
        gnomad_data = None
        

        if variant._genome.name == GenomeReference.hg19.name:
            #gnomad_conn = self.connect_to_database('data/gnomad2/gnomad_db.sqlite3')
            print("Not possible in the current version please use the hg38 version of the variant")
        elif variant._genome.name == GenomeReference.hg38.name:
            

            try:
                with self._gnomad_db as dbhandler:
                    gnomad_cursor = dbhandler._conn.cursor()
                    
                    #gnomad_query = f'SELECT chr as CHR,pos as POS,ref as REF,alt as ALT, AC, AN, AF FROM gnomad_db WHERE chr = {variant._chr} AND pos = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
                    
                    
                    gnomad_query = """
                        SELECT chr as CHR, pos as POS, ref as REF, alt as ALT, AC, AN, AF
                        FROM gnomad_db
                        WHERE chr = ? AND pos = ? AND ref = ? AND alt = ?
                        """
                    
                    params = (variant._chr, variant._pos, variant._ref, variant._alt)
                    
                    gnomad_cursor.execute(gnomad_query, params)

                    gnomad_rows = gnomad_cursor.fetchall()
                    if len(gnomad_rows) > 0:
                        gnomad_data = gnomad_rows[0]
                        gnomad_data = {
                            "CHR": gnomad_data[0],
                            "POS": gnomad_data[1],
                            "REF": gnomad_data[2],
                            "ALT": gnomad_data[3],
                            "AC": gnomad_data[4],
                            "AN": gnomad_data[5],
                            "AF": gnomad_data[6]
                        }
                    else:
                        logging.error("Variant not found")

            except Error as e:
                logging.error(f"Connection to Gnomad failed")
    
        return gnomad_data

    def get_gnomad_constraints(self,gene=''):
        
        """
        
        Get the constraints from gnomad (synonymous z score, missense z score, loss of function z scores, probability of loss of function intolerance) of a given gene
        
        Args:
            gene: A string representing the HGNC Symbol of a gene
        
        Returns:
            gnomad_data: A dictionary of the synonymous z score, missense z score, loss of function z scores, probability of loss of function intolerance)
        
        Examples:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> gnomad = symetrics.get_gnomad_constraints(gene = 'A1BG')
        
        """

        gnomad_data =  None

        try:
            with self._constraints as dbhandler:
                            
                cursor = dbhandler._conn.cursor()
                            
                query = f"SELECT * FROM GNOMADv4Constraints WHERE gene = ?"
                            
                params = (gene,)

                cursor.execute(query,params)
                rows = cursor.fetchall()
                gnomad_data = rows[0]
                gnomad_data = {
                                "gene": gnomad_data[0],
                                "transcript": gnomad_data[1],
                                "syn_z": gnomad_data[2],
                                "mis_z": gnomad_data[3],
                                "lof_z": gnomad_data[4],
                                "pLI": gnomad_data[5]
                }

        except Error as e:
                    logging.error(e)
                    logging.error(f"Connection to {self._constraints} failed")


    
        return gnomad_data

    def liftover(self,variant: VariantObject):

        """
        
        Perform a conversion of the variant position based from their original reference to a target reference. If hg38 is given, it will
        be converted to hg19 and otherwise

        Args:
            variant: A VariantObject instance representing the chromosome, position, reference allele and alternative allele of a variant
        
        Returns:
            liftover_variant: A VariantObject instance representing the chromosome, position, reference allele and alternative allele of a variant after liftover
        
        Exampless:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> variant_hg19 = VariantObject(chr='7',pos='91763673',ref='C',alt='A',genome=GenomeReference.hg19)
            >>> variant_hg38 = symetrics.liftover(variant_hg19)

        """


        liftover_variant = None

        try:
            # synvep is hg38 (pos_GRCh38) abd hg19 (pos)
            with self._db as dbhandler:
                synvep_cursor = dbhandler._conn.cursor()
                synvep_query = ''
                if variant._genome == GenomeReference.hg38:
                    new_reference = GenomeReference.hg19
                    synvep_query = """
                        SELECT chr as CHR, pos as POS, ref as REF, alt as ALT, HGNC_gene_symbol as GENE, synVep as SYNVEP
                        FROM SYNVEP
                        WHERE chr = ? AND pos_GRCh38 = ? AND ref = ? AND alt = ?
                        """    

                elif variant._genome == GenomeReference.hg19:
                    new_reference = GenomeReference.hg38
                    synvep_query = f'SELECT chr as CHR,pos_GRCh38 as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE chr = {variant._chr} AND pos = {variant._pos} AND ref = "{variant._ref}" AND alt = "{variant._alt}"'
                    synvep_query = """
                        SELECT chr as CHR, pos_GRCh38 as POS, ref as REF, alt as ALT, HGNC_gene_symbol as GENE, synVep as SYNVEP
                        FROM SYNVEP
                        WHERE chr = ? AND pos = ? AND ref = ? AND alt = ?
                        """    
                params = (variant._chr, variant._pos, variant._ref, variant._alt)
                
                synvep_cursor.execute(synvep_query,params)
                synvep_rows = synvep_cursor.fetchall()
                variant_info = synvep_rows[0]
                liftover_variant = VariantObject(
                    chr=variant_info[0],
                    pos=variant_info[1],
                    ref=variant_info[2],
                    alt=variant_info[3],
                    genome=new_reference
                )
            
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return liftover_variant
    

    def get_variant_list(self,gene: str):

        """
        
        Get a list of variants associated with the gene

        Args:
            gene: A string type input referring to the gene of interest
        
        Returns:
            variant_list: A list of all variants under the gene
                
        Exampless:

            >>> from symetrics import *
            >>> symetrics = Symetrics('symetrics.db')
            >>> variant_list = symetrics.get_variant_list('A1BG')

        """


        variant_list = None

        try:
            # synvep is hg38 (pos_GRCh38) abd hg19 (pos)
            with self._db as dbhandler:
                synvep_cursor = dbhandler._conn.cursor()
                synvep_query = ''
                synvep_query = f'SELECT chr as CHR,pos as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE HGNC_gene_symbol = "{gene}"'                    
                
                query = "SELECT chr as CHR,pos as POS,ref as REF,alt as ALT, HGNC_gene_symbol as GENE,synVep as SYNVEP FROM SYNVEP WHERE HGNC_gene_symbol = ?"
                            
                params = (gene,)
                
                
                synvep_cursor.execute(query,params)
                synvep_rows = synvep_cursor.fetchall()
                variant_list = synvep_rows
            
        except Error as e:
            logging.error(f"Connection to {self._db} failed")
    
        return variant_list
        

