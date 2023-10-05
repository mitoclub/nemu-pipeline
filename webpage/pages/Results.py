import streamlit as st

st.set_page_config(
    page_title="NeMu Results",
    page_icon="🦈",
    layout="centered",
    initial_sidebar_state="collapsed",
    menu_items={
        'About': "# This is a header. This is an *extremely* cool app!",
        'Get Help': 'https://github.com/mitoclub/nemu-pipeline/wiki',
        'Report a bug': "https://github.com/mitoclub/nemu-pipeline/issues",
    }
)

st.title("NeMu results")
st.sidebar.header("NeMu Results")

params = st.experimental_get_query_params()
st.text('experimental_get_query_params:')
st.json(params)
for key in ['job_id', 'email']:
    if key not in st.session_state and key in params:
        st.session_state[key] = params[key][0]

st.subheader('st.session_state:')
st.json(st.session_state)
