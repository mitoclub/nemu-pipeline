_GENETIC_CODES = [1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33]
_GENETIC_CODES_MITO = [2,3,4,5,9,13,14,16,21,22,23,24,33]
DBS = ['CO1','CO2','CO3','Cytb','A6','A8','ND1','ND2','ND3','ND4','ND4L','ND5','ND6']
SAMPLING_AUTO, SAMPLING_CUSTOM = 'Auto','Custom'
COMPARATIVE, SPECIES = 'Comparative-species', 'Species-specific'

gencodes = {
    1: '1. The Standard Code',
    2: '2. The Vertebrate Mitochondrial Code',
    3: '3. The Yeast Mitochondrial Code',
    4: '4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code',
    5: '5. The Invertebrate Mitochondrial Code',
    6: '6. The Ciliate, Dasycladacean and Hexamita Nuclear Code',
    9: '9. The Echinoderm and Flatworm Mitochondrial Code',
    10: '10. The Euplotid Nuclear Code',
    11: '11. The Bacterial, Archaeal and Plant Plastid Code',
    12: '12. The Alternative Yeast Nuclear Code',
    13: '13. The Ascidian Mitochondrial Code',
    14: '14. The Alternative Flatworm Mitochondrial Code',
    15: '15. Blepharisma Nuclear Code',
    16: '16. Chlorophycean Mitochondrial Code',
    21: '21. Trematode Mitochondrial Code',
    22: '22. Scenedesmus obliquus Mitochondrial Code',
    23: '23. Thraustochytrium Mitochondrial Code',
    24: '24. Rhabdopleuridae Mitochondrial Code',
    25: '25. Candidate Division SR1 and Gracilibacteria Code',
    26: '26. Pachysolen tannophilus Nuclear Code',
    27: '27. Karyorelict Nuclear Code',
    28: '28. Condylostoma Nuclear Code',
    29: '29. Mesodinium Nuclear Code',
    30: '30. Peritrich Nuclear Code',
    31: '31. Blastocrithidia Nuclear Code',
    33: '33. Cephalodiscidae Mitochondrial UAA-Tyr Code',
}
GENETIC_CODES = list(gencodes.keys())
GENETIC_CODES_MITO = [k for k,v in gencodes.items() if "Mito" in v]

assert GENETIC_CODES == _GENETIC_CODES
assert GENETIC_CODES_MITO == _GENETIC_CODES_MITO


def gencode_id2title(s: str):
    return gencodes.get(s, None)

def is_protein(s):
    return

def contain_outgroup():
    return

def is_nucl_multi_fasta():
    return

def run_pipeline(params):
    ok = True
    msg = ''
    if params['sampling'] == SAMPLING_AUTO:
        for prm, prm_name in [('species_name', 'species name'), ('seqs', 'protein sequence')]:
            if not params[prm]:
                ok, msg = False, f'Please specify {prm_name} in the field above'
                break

        if ok and not is_protein(params['seqs']):
            ok, msg = False, 'Query sequence is not a protein'

    if params['sampling'] == SAMPLING_CUSTOM:
        for prm, prm_name in [('outgrp', 'outgroup name or id'), ('seqs', 'nucleotide sequences fasta file')]:
            if not params[prm]:
                ok = False
                msg = f'Please specify {prm_name} in the field above'
                break

    return ok, msg


CONFIG_NEMU_HEAD = '''
// Process Config:

singularity.enabled = true

process {
  container = '/home/dolphin/image_pipeline.sif'
  executor = 'local'
  cpus = '1'
  memory = '2 GB'
}

executor {
  $local {
    cpus = '4'
    memory = '20 GB'
  }  
}
'''

CONFIG_NEMU_BODY = '''
// Run Parameters:

params.sequence = 'sequences.fasta'
params.gencode = '{gencode}'
params.nspecies = '{nspecies}'
params.outgroup = '{outgrp}'
params.aligned = '{aligned}'
params.species_name = '{species_name}'
params.Mt_DB = '/db/MIDORI2_UNIQ_NUC_SP_GB253_{gene_db}_BLAST'
params.verbose = 'true'

// Process Parameters:

// Process Parameters for raxml_build_tree:
params.raxml_build_tree.run_RAXML = "{run_raxml}" //* @checkbox @description:"Check the box to use RAXML for phylogeny reconstruction."
params.raxml_build_tree.raxml_model = "{model_raxml}" //* @input @description:"Substitution model for RAxML"

// Process Parameters for iqtree_build_tree:
params.iqtree_build_tree.run_IQTREE = "{run_iqtree}" //* @checkbox @description:"Checkk the box to use RAXML for phylogeny reconstruction."
params.iqtree_build_tree.iqtree_model = "{model_iqtree}" //* @input @description:"Substitution model for IQTREE2. Could be 'MFP' to run modelfinder."
params.iqtree_build_tree.quantile = "{shrink_quantile}" //* @input @description:"The quantile(s) to set threshold for treeshrink."

// Process Parameters for Phylogeny:
params.phylo.run_shrinking = "{run_shrinking}"

// Process Parameters for mut_processing_params:
params.mut_processing_params.syn4f = "{syn4f}" //* @checkbox @description:"Run extraction of mutational spectrum based on synonymous fourfold mutations"
params.mut_processing_params.all = "{mall}" //* @checkbox @description:"Run extraction of mutational spectrum based on all mutations"
params.mut_processing_params.mnum192 = "{mnum192}" //* @input @description:"Number of mutation types (max 192) required to calculate and plot 192-component mutational spectra"
params.mut_processing_params.use_probabilities = "{use_proba}" //* @checkbox @description:"Use probabilities of nucleotides in mutational spectra calculation"
params.mut_processing_params.proba_min = "{proba_cutoff}" //* @input @description:"Mutation probability cutoff: mutations with lower probability will not be considered in spectra calculation"
params.mut_processing_params.run_pyvolve = "{simulate}" //* @checkbox @description:"Run simulation using pyvolve. Available only when nspecies = 'multiple'"
params.mut_processing_params.replics = "{nreplics}" //* @input @description:"Number of replics to simulate neutral evolution in pyvolve"
params.mut_processing_params.scale_tree = "{scale_tree}" //* @input @description:"Scaling coefficient for tree in pyvolve: less"
'''

def prepare_config(map):
    lines = []
    for line in CONFIG_NEMU_BODY.split('\n'):
        new_line = line.format_map(map)
        lines.append(new_line)
    
    return CONFIG_NEMU_HEAD + '\n'.join(lines)


def prepare_sequences(path):
    # TODO
    return
