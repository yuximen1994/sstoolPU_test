import streamlit as st
import pandas as pd
import numpy as np
import math
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt
from ssfunc import pf_calc
from ssfunc import lsm_sys
from ssdata import case_3bus
from itertools import product

st.set_page_config(layout="wide")

vector1 = ['GFM_Droop', '"GFM_VSM', 'GFL', 'SG']
vector2 = ['GFM_Droop', '"GFM_VSM', 'GFL', 'SG']
# Create all combinations and index them from 'com1' to 'com16'
combinations = list(product(vector1, vector2))
combination_named_index = {f'com{index + 1}': combination for index, combination in enumerate(combinations)}
selected_combination1 = combination_named_index['com1']
selected_combination2 = combination_named_index['com2']
# sidebar
sidebar1 = st.sidebar.selectbox(
    "What configuration do you want to select for the 1st generator?",
    ("GFM_Droop", "GFM_VSM", "GFL", "SG"), index=0, placeholder="Select configuration...",
)
sidebar2 = st.sidebar.selectbox(
    "What configuration do you want to select for the 2nd generator?",
    ("GFM_Droop", "GFM_VSM", "GFL", "SG"), index=0, placeholder="Select configuration...",
)
if (sidebar1,sidebar2) == ('GFM_Droop','GFM_Droop'):
    st.write('You selected GFM_Droop and GFM_Droop.')
else:
    st.write("The selected combination is not supported.")

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
numeigs = len(eigvals)
lefteigenvectors = np.linalg.inv(eigenvectors)
pmatrix = np.multiply(eigenvectors,np.transpose(lefteigenvectors))
pmatrixabs = abs(pmatrix)
stateVariableNames = ['theta1','P01','Qo1','phid1','phiq1','gammad1','gammaq1','iid1','iiq1','vcd1','vcq1','iod1','ioq1',
                     'theta2','epsilonL2','wf2','P02','Qo2','phid2','phiq2','gammad2','gammaq2','iid2','iiq2','vcd2','vcq2','iod2','ioq2',
                     'ibranchD1','ibranchQ1','ibranchD2','ibranchQ2','iloadD','iloadQ']
modeNames = ['mode{}'.format(i) for i in range(1,numeigs+1)]
# plot participation factor map
figheatmap = px.imshow(pmatrixabs,
                       labels=dict(x="modes", y="state variables"),
                       x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34],
                       y = stateVariableNames)
figheatmap.update_layout(height=800)

# Use text_input for manual number input
input_number = st.sidebar.text_input("Which mode do you want to select? (1-"+str(numeigs)+")", value='1')
number = int(input_number)
# Check if the input is a number and within the desired range
if input_number:
    try:        
        number = int(input_number) # Convert input to an integer        
        if 1 <= number <= numeigs: # Check if the number is in the range
            df = pd.DataFrame(pmatrixabs, columns=modeNames)
            df.insert(0, "statevariables", stateVariableNames, True)           
            df.loc[df[modeNames[number-1]] < 0.02, 'statevariables'] = 'Other states'  # Represent state variables with a relatively larger participation factor
            figpie = px.pie(df, values=modeNames[number-1], names='statevariables', title='Participation factor analysis of mode '+str(number))
            figpie.update_layout(title={'text':'Participation factor analysis of mode '+str(number),'x':0.415,'xanchor':'center'})
            st.text("real: ")
            st.text("imag: ")
        else:
            st.error('Number out of range. Please enter a number between 1 and '+str(numeigs)+'.')
    except ValueError:        
        st.error('Invalid input. Please enter a number.') # Handle the case where input is not a number

col1, col2 = st.columns(2,gap="small")
with col1:
   st.plotly_chart(figheatmap, height=800, theme="streamlit",use_container_width=True)
with col2:
   st.plotly_chart(figpie, height=800, theme="streamlit",use_container_width=True)

# plot table
mode = range(1,len(eigvals)+1)
realpart = eigvals.real
imagpart = eigvals.imag
frequency = eigvals.imag/2/math.pi
dampingratio = -eigvals.real/np.sqrt(realpart*realpart+imagpart*imagpart)
list_of_tuples = list(zip(mode, realpart, imagpart, frequency, dampingratio)) 
df = pd.DataFrame(list_of_tuples,
                 columns = ["mode", "real", "image", "Frequency(Hz)", "damping ratio"])
st.table(df)
