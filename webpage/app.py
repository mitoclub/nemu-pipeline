'''
usage: streamlit run demo.py

Page scheme:
0. done description 
1. sampling type (auto, custom data)
2. inter- or intra-species mode
3. main params
4. hided optional params with defaults
5. results subpage that will show images and suggest to download results
'''

import random

import streamlit as st

from utils import GENETIC_CODES_MITO, GENETIC_CODES, gencode_id2title


st.set_page_config(
    page_title="NeMu pipeline",
    page_icon="🦈",
    layout="centered",
    # initial_sidebar_state="expanded",
    menu_items={
        'About': "# This is a header. This is an *extremely* cool app!",
        'Get Help': 'https://github.com/mitoclub/nemu-pipeline/wiki',
        'Report a bug': "https://github.com/mitoclub/nemu-pipeline/issues",
    }
)
st.title('NeMu pipeline')
st.markdown('The pipeline for neutral mutation spectra evaluation based on evolutionary data for one species')
st.markdown('Execute pipeline on IKBFU server')


dbs = ['CO1','CO2','CO3','Cytb','A6','A8','ND1','ND2','ND3','ND4','ND4L','ND5','ND6']
SAMPLING_AUTO, SAMPLING_CUSTOM = 'Auto','Custom'
COMPARATIVE, SPECIES = 'Comparative-species', 'Species-specific'

st.header("Inputs")
job_title = st.text_input("Job Title", placeholder='Run-345', max_chars=100)
sampling = st.radio('Sampling type', [SAMPLING_AUTO, SAMPLING_CUSTOM], help=f'Select "{SAMPLING_AUTO}" to search homologous sequences in the blast-dbs of mitochondrial genes or "{SAMPLING_CUSTOM}" to execute pipeline on custom homologous sequences')
if sampling == SAMPLING_CUSTOM:
    level = st.radio('Level', [SPECIES, COMPARATIVE], help=f'"{SPECIES}" to analyse single species, "{COMPARATIVE}" to analyse several species')
    gencode = st.selectbox('Genetic code', GENETIC_CODES, 1, format_func=gencode_id2title)
    gene_db = species_name = None
    seqs = st.file_uploader('File uploader')
    outgrp = st.text_input("Outgroup", placeholder="OUTGRP", help="Id or header of outgroup record in query fasta (>Id)")
    
elif sampling == SAMPLING_AUTO:
    level = st.radio('Level of analysis', [SPECIES], help=f'"{COMPARATIVE}" to analyse several species, "{SPECIES}" to analyse single species')
    gencode = st.selectbox('Genetic code', GENETIC_CODES_MITO, format_func=gencode_id2title)
    gene_db = st.selectbox('mtDNA gene', dbs, help='Select precomputed blast-db from MIDORI2 reference database of mitochondrial DNA sequences based on GenBank version 253 (GB253)')
    species_name = st.text_input('Species name')
    seqs = st.text_area("Query fasta", max_chars=10000, placeholder=">id name\nACGTCCTG...")
    outgrp = None

with st.expander("**Advanced pipeline parameters**"):
    col1, col2 = st.columns(2)
    with col1:
        run_iqtree = st.checkbox('Use IQTREE2', True, help='Check the box to use IQTREE2 for phylogeny reconstruction (http://www.iqtree.org/doc/Substitution-Models)')
        model_iqtree = '10.12+FO+G6+I' if level == COMPARATIVE else 'GTR+FO+G6+I'
        if run_iqtree:
            model_iqtree = st.text_input('Substitution model for IQTREE2', model_iqtree)
        shrink_quantile = st.number_input("Quantile for TreeShrink", 0., 1., 0.05, 0.01)

    with col2:
        run_raxml = st.checkbox('Use RAxML', False, help='Check the box to use RAXML for phylogeny reconstruction')
        model_raxml = None
        if run_raxml:
            model_raxml = st.text_input('Substitution model for RAxML', 'GTRGAMMAIX')

    syn4f = st.checkbox('Synonymous fourfold mutations', help="Run extraction of mutational spectrum based on synonymous fourfold mutations")
    mall = st.checkbox('Synonymous and nonsynonymous', help='Run extraction of mutational spectrum based on all mutations: synonymous and nonsynonymous')
    proba = st.checkbox('Use probabilities', True, help='Use probabilities of nucleotides in mutational spectra calculation')
    col3, col4 = st.columns(2)
    with col3:
        proba_cutoff = st.number_input('Mutation probability cutoff', 0., 1., 0.3, 0.05, help="Mutation probability cutoff: mutations with lower probability will not be considered in spectra calculation")
    with col4:
        mnum192 = st.number_input("Mutation types cutoff", 1, 192, 16, 1, help='Number of mutation types (max 192) required to calculate and plot 192-component mutational spectrum')
    
    if level == SPECIES:
        simulate = st.checkbox('Run simulation of neutral evolution', help='Run simulation using pyvolve to estimate mutation spectra neutrality. Available only for species-specific analysis')
        if simulate:
            col5, col6 = st.columns(2)
            with col5:
                nreplics = st.number_input("Number of replics", 1, 50, 10, 1, help='Number of replics to simulate neutral evolution in pyvolve')
            with col6:
                scale_tree = st.number_input("Tree scaling coefficient", 0., 10., 1., help='Scaling coefficient for tree in pyvolve')
        else:
            nreplics = scale_tree = None

    else:
        simulate = False
    site_rates = st.checkbox('Run site rates estimation in phylogenetic inference')
    if site_rates:
        col7, col8 = st.columns(2)
        with col7:
            ncategories = st.number_input("Number of categories", 2, 8, 6, 1, help='Site rates from discrete Gamma model with N categories')
        with col8:
            cat_cutoff = st.number_input("Category cutoff", 0, ncategories, 1, 1, help='Minimal category of sites that will be used in synonymous nucleoride frequencies calculation')
    else:
        ncategories = cat_cutoff = None

params = {
    'job_title': job_title if job_title else f'Run-{random.randint(100, 10000)}',
    'sampling': sampling,
    'level': level,
    'species_name': species_name,
    'gencode': gencode,
    'gene_db': gene_db,
    'species_name': species_name,
    'seqs': seqs,
    'outgrp': outgrp,
    'run_iqtree': run_iqtree,
    'model_iqtree': model_iqtree,
    'run_raxml': run_raxml,
    'model_raxml': model_raxml,
    'shrink_quantile': shrink_quantile,
    'syn4f': syn4f,
    'mall': mall,
    'proba': proba,
    'proba_cutoff': proba_cutoff,
    'mnum192': mnum192,
    'simulate': simulate,
    'nreplics': nreplics,
    'scale_tree': scale_tree,
    'site_rates': site_rates,
    'ncategories': ncategories,
    'cat_cutoff': cat_cutoff,


}
st.json(params)
st.button('Run pipeline!')


# with st.sidebar:
#     with st.echo():
#         st.write("This code will be printed to the sidebar.")

#     with st.spinner("Loading..."):
#         time.sleep(5)
#     st.success("Done!")


# on = st.toggle('Activate feature')
# if on:
#     st.write('Feature activated!')
