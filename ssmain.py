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

fig = px.imshow(abs(pmatrix),
                labels=dict(x="modes", y="state variables"),
                x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]
                y = ['theta1','P01','Qo1','phid1','phiq1','gammad1','gammaq1','iid1','iiq1','vcd1','vcq1','iod1','ioq1',
                     'theta2','epsilonL2','wf2','P02','Qo2','phid2','phiq2','gammad2','gammaq2','iid2','iiq2','vcd2','vcq2','iod2','ioq2',
                     'ibranchD1','ibranchQ1','ibranchD2','ibranchQ2','iloadD','iloadQ']
                )
fig.update_layout(height=800)
st.plotly_chart(fig, height=800, theme="streamlit")
