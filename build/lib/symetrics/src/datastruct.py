from enum import Enum,auto

class MetricsGroup(Enum):
    '''
    Lists of all possible metrics for the package
    
    '''
    SURF = auto()
    RSCU = auto()
    CPGX = auto()
    CPG = auto()
    DRSCU = auto()
    GERP = auto()
    SYNVEP = auto()
    SPLICEAI = auto()

    SPLICE = auto()
    SYMETRICS = auto()
    MES = auto()
    CPGEXON = auto()
    CPGLogit = auto()

class GenomeReference(Enum):

    '''
    Lists of all possible genome reference for the package
    
    '''
    hg19 = auto()
    hg38 = auto()
    
class VariantObject():

    '''
    An object representing the chromosome, position, reference allele and alternative allele of a variant 
    
    '''

    _chr = ''
    _pos = ''
    _alt = ''
    _ref = ''
    _genome = 'hg38'

    def __init__(self,chr:str = '',pos:str = '', ref:str = '',alt:str = '',genome:GenomeReference = GenomeReference.hg38) -> None:
        self._chr = chr
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._genome = genome
        