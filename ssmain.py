import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from ssfunc import pf_calc
from ssfunc import lsm_sys
from ssdata import case_3bus

chart_data = pd.DataFrame(np.random.randn(20, 3), columns=["a", "b", "c"])

######################################
# Define your images
images = {
     'GFM_Droop-GFM_Droop': 'fig/DroopDroop.png',
     'GFM_Droop-GFM_VSM': 'fig/DroopVSM.png',
     'GFM_Droop-GFL': 'fig/DroopGFL.png',
     'GFM_Droop-SG': 'fig/SGDroop.png',
     'GFM_VSM-GFM_Droop': 'fig/DroopVSM.png',
     'GFM_VSM-GFM_VSM': 'fig/VSMVSM.png',
     'GFM_VSM-GFL': 'fig/VSMGFL.png',
     'GFM_VSM-SG': 'fig/SGVSM.png',
     'GFL-GFM_Droop': 'fig/DroopGFL.png',
     'GFL-GFM_VSM': 'fig/VSMGFL.png',
     'GFL-SG': 'fig/SGGFL.png',
     'SG-GFM_Droop': 'fig/SGDroop.png',
     'SG-GFM_VSM': 'fig/SGVSM.png',
     'SG-GFL': 'fig/SGGFL.png',
     'SG-SG': 'fig/SGSG.png',
}

# Create the selectboxes
selectbox1 = st.selectbox('What control strategy do you want to select for the 1st inverter?', ['GFM_Droop', 'GFM_VSM', 'GFL', 'SG'], index=None, placeholder="Select control method...",)
selectbox2 = st.selectbox('What control strategy do you want to select for the 2nd inverter?', ['GFM_Droop', 'GFM_VSM', 'GFL', 'SG'], index=None, placeholder="Select control method...",)

# Determine which image to display
selected_image = images.get(f'{selectbox1}-{selectbox2}')

# Display the image
if selected_image:
    st.image(selected_image)
else:
    st.write("No control diagram available.")
######################################


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
                x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34],
                y = ['theta1','P01','Qo1','phid1','phiq1','gammad1','gammaq1','iid1','iiq1','vcd1','vcq1','iod1','ioq1',
                     'theta2','epsilonL2','wf2','P02','Qo2','phid2','phiq2','gammad2','gammaq2','iid2','iiq2','vcd2','vcq2','iod2','ioq2',
                     'ibranchD1','ibranchQ1','ibranchD2','ibranchQ2','iloadD','iloadQ']
                )
fig.update_layout(height=800)
st.plotly_chart(fig, height=800, theme="streamlit")

mode = range(1,len(eigvals))
realpart = eigvals.real
imagpart = eigvals.imag
frequency = eigvals.imag/2/math.pi
dampingratio = -eigvals.real/eigvals.real
list_of_tuples = list(zip(mode, realpart, imagpart, frequency, dampingratio)) 
df = pd.DataFrame(list_of_tuples,
                 columns = ["mode", "real", "image", "Frequency(Hz)", "damping ratio"])


st.table(df)
