import streamlit as st
import pandas as pd
import numpy as np
from ssfunc import pf_calc
from ssfunc import lsm_sys
from ssdata import case_3bus

chart_data = pd.DataFrame(np.random.randn(20, 3), columns=["a", "b", "c"])

st.bar_chart(chart_data)
