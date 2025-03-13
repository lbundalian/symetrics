from symetrics import *
import os
import time

cfg = os.path.join(os.path.dirname(__file__), 'config.json')
api = Symetrics(cfg)

chr_val = 'X'
pos_val = '153906421 '
ref_val = 'G'
alt_val = 'C'
assy_val = 'hg38'

variant = VariantObject(chr=chr_val,pos=pos_val,ref=ref_val,alt=alt_val,genome=GenomeReference.hg38)
assy = "hg38"


t = api.get_prop_score(MetricsGroup.SPLICEAI,"GPR39")

s = time.time()

### start
liftover_variant = api.liftover(variant)

if assy == 'hg38':
        variant_hg19, variant_val = liftover_variant, variant
else:
        variant_hg19, variant_val = variant, liftover_variant

# Fetch metric scores
metrics = {
        'Synvep': api.get_synvep_score(variant_val),
        'SPLICEAI': api.get_spliceai_score(variant=variant_val),
        'SURF': api.get_surf_score(variant=variant_val),
        'SILVA': api.get_silva_score(variant=variant_hg19),
        'GNOMAD': api.get_gnomad_data(variant=variant_val)
}

# t = api.get_prop_score(MetricsGroup.SYNVEP.name,'A1BG')
# Helper function for safe float conversion
def safe_float(value):
        return float(value) if value and value != 0 else 0

# Prepare model input
silva = metrics['SILVA']
model_input = {
        'Synvep': safe_float(metrics['Synvep'].get('SYNVEP')),
        'SPLICEAI': safe_float(metrics['SPLICEAI'].get('MAX_DS')) if metrics['SPLICEAI'] else 'N/A',
        'SURF': safe_float(metrics['SURF'].get('SURF')),
        'MES': safe_float(silva.get('MES')),
        'GERP': safe_float(silva.get('GERP')),
        'CpG': safe_float(silva.get('CPG')),
        'CpG_exon': safe_float(silva.get('CPGX')),
        'RSCU': safe_float(silva.get('RSCU')),
        'dRSCU': safe_float(silva.get('dRSCU')),
        'F_MRNA': safe_float(silva.get('F_MRNA')),
        'F_PREMRNA': safe_float(silva.get('F_PREMRNA')),
        'AF': metrics['GNOMAD'].get('AF', 0)
}

print(model_input)
# Predict probability
probability = api.predict_probability(model_input)
deleterious_prob = round(probability[0][1], 2) 
# Prepare result output
results = {
        'Information': {
        'Gene': metrics['Synvep'].get('GENE', 'N/A'),
        'Variants': f"chr{chr_val}:{pos_val} {ref_val}>{alt_val}"
        },
        'Metrics': {key: {'value': val, 'threshold': threshold, 'description': description}
                for key, val, threshold, description in [
                        ('Synvep', model_input['Synvep'], 0.5, "Probability of human synonymous variant to have an effect/no effect. (t = 0.5)"),
                        ('SPLICEAI', model_input['SPLICEAI'], 0.5, "Uses the Max Delta Scores which is the maximum possible delta scores quantifying how much variant is expected to alter splicing pattern. (t = 0.5)"),
                        ('SURF', model_input['SURF'], 5, "Measure of RNA-folding and stability to capture the effects of synonymous variation to secondary mRNA structure. (t = 5)"),
                        ('MES', model_input['MES'], 3, "Maximum efficiency of splicing site motif or junction splice strength. (t = 3)"),
                        ('GERP', model_input['GERP'], 4, "Measure of conservation for a given synonymous mutation by calculating site-specific rejected substitution scores. (t = 4)"),
                        ('CpG', model_input['CpG'], 0, "Observed/expected CpG content of exon (t = 0)"),
                        ('CpG_exon', model_input['CpG_exon'], 1, "Binary values implying whether a mutation changes a CpG or not (t = 1)"),
                        ('RSCU', model_input['RSCU'], 1, "Relative synonymous codon usage of new codon (t = 1)"),
                        ('dRSCU', model_input['dRSCU'], 1, "Change in Relative synonymous codon usage caused by mutation (t = 1)"),
                        ('F_MRNA', model_input['F_MRNA'], 1, "Relative distance to the end of mature MRNA (t = 1)"),
                        ('PREMRNA', model_input['F_PREMRNA'], 1, "Relative distance to the end of pre-MRNA (t = 1)")
                ]},
        'Gnomad': {
        'AC': metrics['GNOMAD'].get('AC', 0),
        'AN': metrics['GNOMAD'].get('AN', 0),
        'AF': metrics['GNOMAD'].get('AF', 0)
        },
        'Probability': {'DELETERIOUS': deleterious_prob}
}
### end


e = time.time()
print(e - s)

