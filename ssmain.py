import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from ssfunc import pf_calc
from ssfunc import lsm_sys
from ssdata import case_3bus

chart_data = pd.DataFrame(np.random.randn(20, 3), columns=["a", "b", "c"])

sysData = case_3bus()
x, xdot = lsm_sys(sysData)
Xss = [0.0000,0.5147,0.1411,0.1452,-0.0386,0.0844,-0.0002,0.5228,
       -0.1388,0.0386,-0.1930,0.5132,-0.1407,0.0000,0.0000,0.9993,
       0.5147,0.1411,0.1452,-0.0386,0.0844,-0.0002,0.5228,-0.1388,
       0.0386,-0.1930,0.5132,-0.1407,1.0249,-0.2814,0.5127,-0.1407,
       0.5127,-0.1407]
Asys = xdot.jacobian(x)
for i in range(len(Xss)):
    Asys = Asys.subs([(x[i], Xss[i])])

Asys = np.array(Asys).astype(np.float64)
eigvals, eigenvectors = np.linalg.eig(Asys)
lefteigenvectors = np.linalg.inv(eigenvectors)
pmatrix = np.multiply(eigenvectors,np.transpose(lefteigenvectors))

fig = px.imshow(abs(pmatrix))
st.plotly_chart(fig, theme="streamlit")
