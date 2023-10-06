'''
usage: streamlit run demo.py

Page scheme:
0. done description 
1. done sampling type (auto, custom data)
2. done inter- or intra-species mode
3. done main params
4. done hided optional params with defaults
5. results subpage that will show images and suggest to download results
'''

import random

from dotenv import load_dotenv
import streamlit as st
# import streamlit.config as sconfig

from utils import (
    GENETIC_CODES_MITO, GENETIC_CODES, DBS, 
    SAMPLING_AUTO, SAMPLING_CUSTOM, COMPARATIVE, SPECIES,
    gencode_id2title, run_pipeline, prepare_config,
)
from connection import send_email

load_dotenv()

if 'job_id' not in st.session_state:
    st.session_state.job_id = random.randint(1000, 1000000)

st.set_page_config(
    page_title="NeMu pipeline",
    page_icon="🦈",
    layout="centered",
    initial_sidebar_state="collapsed",
    menu_items={
        'About': "# NeMu pipeline interface",
        'Get Help': 'https://github.com/mitoclub/nemu-pipeline/wiki',
        'Report a bug': "https://github.com/mitoclub/nemu-pipeline/issues",
    }
)
st.title('NeMu pipeline')
st.sidebar.header("NeMu Pipeline")
st.markdown('The pipeline for neutral mutation spectra evaluation based on evolutionary data for one species')
st.markdown('Execute pipeline on IKBFU server')


st.header("Inputs")
job_title = st.text_input("Job Title", placeholder='Test Run (optional)', max_chars=100)
sampling = st.radio('Sampling type*', [SAMPLING_AUTO, SAMPLING_CUSTOM], help=f'Select "{SAMPLING_AUTO}" to search homologous sequences in the blast-dbs of mitochondrial genes or "{SAMPLING_CUSTOM}" to execute pipeline on custom homologous sequences')
if sampling == SAMPLING_CUSTOM:
    level = st.radio('Level of analysis*', [SPECIES, COMPARATIVE], help=f'"{SPECIES}" to analyse single species, "{COMPARATIVE}" to analyse several species')
    gencode = st.selectbox('Genetic code*', GENETIC_CODES, 1, format_func=gencode_id2title, help='Genetic code for alighning and mutations annotation')
    gene_db = species_name = None
    # st.divider()
    outgrp = st.text_input("Outgroup*", placeholder="OUTGRP", help="Id or header of outgroup record in query fasta (>OUTGRP)")
    seqs = st.file_uploader('Query fasta*', ['fasta', 'fna', 'fa'], help='Fasta-file with several nucleotide sequences')

elif sampling == SAMPLING_AUTO:
    level = st.radio('Level of analysis*', [SPECIES], help=f'"{COMPARATIVE}" to analyse several species, "{SPECIES}" to analyse single species')
    gencode = st.selectbox('Genetic code*', GENETIC_CODES_MITO, format_func=gencode_id2title, help='Genetic code for blasting, alighning and mutations annotation')
    # st.divider()
    # https://docs.streamlit.io/library/advanced-features/forms is shit
    col1, col2 = st.columns(2)
    with col1:
        species_name = st.text_input('Species name*',  max_chars=100, help='For example *Canis lupus*')
    with col2:
        gene_db = st.selectbox('mtDNA gene*', DBS, help='Select precomputed blast-db from MIDORI2 reference database of mitochondrial DNA sequences based on GenBank version 253 (GB253)')
    seqs = st.text_area("Query protein sequence*", max_chars=5000, placeholder="TSKHHFGFQAAAWYWHFVDVVWLFLYVSIYWWGS...", help='Single protein sequence that will be used in homoloug nucleotide sequences search using tblastn ')
    outgrp = None

with st.expander("**Advanced pipeline parameters**"):
    # TODO add sections: phylo params ...
    col1, col2 = st.columns(2)
    with col1:
        run_iqtree = st.checkbox('Use IQTREE2', True, help='Check the box to use IQTREE2 for phylogeny reconstruction (http://www.iqtree.org/doc/Substitution-Models)')
        model_iqtree = '10.12+FO+G6+I' if level == COMPARATIVE else 'GTR+FO+G6+I'
        if run_iqtree:
            model_iqtree = st.text_input('Substitution model for IQTREE2', model_iqtree)
        run_shrinking = st.checkbox('Use Tree Shrinking', True, help='Filter out outlier nodes (https://uym2.github.io/TreeShrink/)')
        if run_shrinking:
            shrink_quantile = st.number_input("Quantile for TreeShrink", 0., 1., 0.05, 0.01)

    with col2:
        run_raxml = st.checkbox('Use RAxML', False, help='Check the box to use RAXML for phylogeny reconstruction')
        model_raxml = None
        if run_raxml:
            model_raxml = st.text_input('Substitution model for RAxML', 'GTRGAMMAIX')

    syn4f = st.checkbox('Synonymous fourfold mutations', help="Run extraction of mutational spectrum based on synonymous fourfold mutations")
    mall  = st.checkbox('Synonymous and nonsynonymous', help='Run extraction of mutational spectrum based on all mutations: synonymous and nonsynonymous')
    use_proba = st.checkbox('Use probabilities', True, help='Use probabilities of nucleotides in mutational spectra calculation')
    col3, col4 = st.columns(2)
    with col3:
        if use_proba:
            proba_cutoff = st.number_input('Mutation probability cutoff', 0., 1., 0.3, 0.05, help="Mutation probability cutoff: mutations with lower probability will not be considered in spectra calculation")
        mnum192 = st.number_input("Mutation types cutoff", 1, 192, 16, 1, help='Number of mutation types (max 192) required to calculate and plot 192-component mutational spectrum')

    nreplics = scale_tree = ncategories = cat_cutoff = None
    simulate = False
    if level == SPECIES:
        simulate = st.checkbox('Run simulation of neutral evolution', help='Run simulation using pyvolve to estimate mutation spectra neutrality. Available only for species-specific analysis')
        if simulate:
            col5, col6 = st.columns(2)
            with col5:
                nreplics = st.number_input("Number of replics", 1, 50, 10, 1, help='Number of replics to simulate neutral evolution in pyvolve')
            with col6:
                scale_tree = st.number_input("Tree scaling coefficient", 0., 10., 1., help='Scaling coefficient for tree in pyvolve')

    site_rates = st.checkbox('Run site rates estimation in phylogenetic inference')
    if site_rates:
        col7, col8 = st.columns(2)
        with col7:
            ncategories = st.number_input("Number of categories", 2, 8, 6, 1, help='Site rates from discrete Gamma model with N categories')
        with col8:
            cat_cutoff = st.number_input("Category cutoff", 0, ncategories, 1, 1, help='Minimal category of sites that will be used in synonymous nucleoride frequencies calculation')

email = st.text_input("email", key='email', placeholder='test@example.com', help='Pipeline output will be sent to the specified email')

if email:
    st.button("Send test email", on_click=send_email, 
            args=(st.session_state.email, ['requirements.txt', '.streamlit/config.toml']))

params = {
    'job_id': st.session_state.job_id,
    'job_title': job_title,
    'email': email,
    'sampling': sampling,
    'level': level,
    'nspecies': 'single' if level == SPECIES else 'multiple',
    'aligned': 'TODO',
    'species_name': species_name,
    'gencode': gencode,
    'gene_db': gene_db,
    'species_name': species_name,
    'seqs': seqs,
    'outgrp': outgrp,
    'run_iqtree': str(run_iqtree).lower(),
    'model_iqtree': model_iqtree,
    'run_raxml': str(run_raxml).lower(),
    'model_raxml': model_raxml,
    'run_shrinking': str(run_shrinking).lower(),
    'shrink_quantile': shrink_quantile,
    'syn4f': str(syn4f).lower(),
    'mall': str(mall).lower(),
    'use_proba': str(use_proba).lower(),
    'proba_cutoff': proba_cutoff,
    'mnum192': mnum192,
    'simulate': str(simulate).lower(),
    'nreplics': nreplics,
    'scale_tree': str(scale_tree).lower(),
    'site_rates': site_rates,
    'ncategories': ncategories,
    'cat_cutoff': cat_cutoff,

}

def click_pipeline_btn():
    ok, msg = run_pipeline(params)
    if ok:
        st.text(prepare_config(params))
        st.link_button("Go to Results page", 
                    f"/Results/?job_id={st.session_state.job_id}&email={st.session_state.email}")
    else:
        st.toast(msg, icon='🔥')

    # TODO почему-то все идет вверх страницы, а не на текущем уровне

run_it = st.button('Run pipeline!', on_click=click_pipeline_btn)

# https://github.com/blackary/st_pages ????

# st.markdown('<a href="/Results" target=_top>Results</a>', unsafe_allow_html=True)

st.subheader('st.session_state:')
st.json(st.session_state)

st.subheader('Pipeline params:')
st.json(params)

st.divider()

#### !!!!!!!!!!!!!!
# job_id send to url params and areate new link that will follow to the results page that will read params and show results of current work
# https://docs.streamlit.io/library/api-reference/utilities/st.experimental_get_query_params




# with st.sidebar:
#     with st.echo():
#         st.write("This code will be printed to the sidebar.")

#     with st.spinner("Loading..."):
#         time.sleep(5)
#     st.success("Done!")


# on = st.toggle('Activate feature')
# if on:
#     st.write('Feature activated!')
