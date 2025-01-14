from setuptools import setup, find_packages

VERSION = '15.0.0' 
DESCRIPTION = 'Symetrics API'
LONG_DESCRIPTION = '''A package to access SYMETRICS database. 
SYMETRICS database is a consolidation of metrics for synonymous variants which were derived from a number of computational tools each of which contributing to 
attribute specific metrics such as SYNVEP for general functional constraints, SpliceAI for splicing effect, SILVA for obtaining GERP++ (phylogenetic related constraints)
CpG/CpG_Exon, dRSCU/RSCU for codon usage and SURF for rna stability. The package also includes a result of the analysis of the influence of each variants exceeding set threshold
per metrics defined constituting to a score assigned to a gene'''
REQUIRED_PACKAGES = [
    'numpy==1.24.4',
    'pandas==2.0.2',
    'scikit-learn==1.2.2'
]

setup(
        name="symetrics", 
        version=VERSION,
        author="Linnaeus Bundalian",
        author_email="linnaeusbundalian@gmail.com",
        description=DESCRIPTION,
        url="https://lbundalian.github.io/symetrics/",
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=REQUIRED_PACKAGES,
        keywords=['python', 'synonymous variants'],
        include_package_data=True,
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)