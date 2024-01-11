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

vector1 = ['GFM_Droop', '"GFM_VSM', 'GFL', 'SG']
vector2 = ['GFM_Droop', '"GFM_VSM', 'GFL', 'SG']

# Create all combinations and index them from 'com1' to 'com16'
combinations = list(product(vector1, vector2))
combination_named_index = {f'com{index + 1}': combination for index, combination in enumerate(combinations)}

selected_combination1 = combination_named_index['com1']
selected_combination2 = combination_named_index['com2']

if combination_named_index == 'selected_combination1':
    st.write('You selected GFM_Droop and GFM_Droop.')
else:
    st.write("You didn\'t select GFM_Droop and GFM_Droop.")

# sidebar
sidebar1 = st.sidebar.selectbox(
    "What configuration do you want to select for the 1st generator?",
    ("GFM_Droop", "GFM_VSM", "GFL", "SG"), index=None, placeholder="Select configuration...",
)
sidebar2 = st.sidebar.selectbox(
    "What configuration do you want to select for the 2nd generator?",
    ("GFM_Droop", "GFM_VSM", "GFL", "SG"), index=None, placeholder="Select configuration...",
)

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
pmatrixabs = abs(pmatrix)
stateVariableNames = ['theta1','P01','Qo1','phid1','phiq1','gammad1','gammaq1','iid1','iiq1','vcd1','vcq1','iod1','ioq1',
                     'theta2','epsilonL2','wf2','P02','Qo2','phid2','phiq2','gammad2','gammaq2','iid2','iiq2','vcd2','vcq2','iod2','ioq2',
                     'ibranchD1','ibranchQ1','ibranchD2','ibranchQ2','iloadD','iloadQ']

fig = px.imshow(pmatrixabs,
                labels=dict(x="modes", y="state variables"),
                x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34],
                y = stateVariableNames)
fig.update_layout(height=800)
st.plotly_chart(fig, height=800, theme="streamlit")

NumElement = len(eigvals)

# Use text_input for manual number input
#input_number = st.sidebar.text_input("Which mode do you want to select? (1-"+str(NumElement)+")")
#number = int(input_number)

labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
sizes = [15, 30, 45, 10]
explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
figpie, axpie = plt.subplots()
axpie.pie(pmatrixabs[:,10], explode=np.zeros(NumElement), labels=stateVariableNames, autopct='%1.1f%%',
          shadow=False, startangle=90)
axpie.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
st.pyplot(figpie)


df = pd.DataFrame(pmatrixabs, columns=['col1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34'])
df.insert(0, "statevariables", stateVariableNames, True)
# Represent state variables with a relatively larger participation factor
#df.loc[df['col1'] >= 0.1, 'statevariables'] = 'Other state variables'
fig = px.pie(df, values='col1', names='statevariables', title='Population of European continent')
st.plotly_chart(fig, height=800, theme="streamlit")

# Check if the input is a number and within the desired range
if input_number:
    try:
        # Convert input to an integer
        number = int(input_number)
        # Check if the number is in the range
        if 1 <= number <= NumElement:
            #fig = px.pie(pmatrix[number,:], values='pop', names='country', title='Population of European continent')
            #st.plotly_chart(fig, use_container_width=True)
            st.write('You entered:', number)
        else:
            st.error('Number out of range. Please enter a number between 1 and '+str(NumElement)+'.')
    except ValueError:
        # Handle the case where input is not a number
        st.error('Invalid input. Please enter a number.')

mode = range(1,len(eigvals)+1)
realpart = eigvals.real
imagpart = eigvals.imag
frequency = eigvals.imag/2/math.pi
dampingratio = -eigvals.real/np.sqrt(realpart*realpart+imagpart*imagpart)
list_of_tuples = list(zip(mode, realpart, imagpart, frequency, dampingratio)) 
df = pd.DataFrame(list_of_tuples,
                 columns = ["mode", "real", "image", "Frequency(Hz)", "damping ratio"])


st.table(df)
