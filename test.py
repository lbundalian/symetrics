from symetrics import *
import os

cfg = os.path.join(os.path.dirname(__file__), 'config.json')
api = Symetrics(cfg)

# 7-91763673-G-A

#n = api.get_variant_list('A1BG')
variant = VariantObject(chr='19',pos=' 10305513',ref='G',alt='T',genome=GenomeReference.hg19)
n = api.get_silva_score(variant)
print(n)

variant = VariantObject(chr='19',pos='10194837',ref='G',alt='T',genome=GenomeReference.hg38)
n = api.get_synvep_score(variant)
print(n)

n = api.get_spliceai_score(variant)
print(n)

n = api.get_surf_score(variant)
print(n)

n = api.get_prop_score(group='SYNVEP',gene='DNMT1')
print(n)

n = api.get_prop_score(group='SURF',gene='DNMT1')
print(n)
