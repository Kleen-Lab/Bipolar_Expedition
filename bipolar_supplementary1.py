#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import math
import pandas as pd
from pathlib import Path


# In[2]:


project_dir = Path.cwd()
data_file = project_dir / 'data' / 'tableoutputfile.xlsx'
df = pd.read_excel(data_file)
df


# In[3]:


freq_order = ['delta', 'theta', 'alpha', 'beta', 'gamma', 'high_gamma']
df['freq'] = pd.Categorical(df['freq'], categories=freq_order, ordered=True)
grouped_df = df.groupby(['freq', 'patient', 'bpd']).agg({'power': 'mean'}).reset_index()
pivot_df = grouped_df.pivot_table(index=['freq', 'patient'], columns=['bpd'], values='power', aggfunc='first')
array_3d = pivot_df.to_numpy().reshape(len(freq_order), len(df['patient'].unique()), len(df['bpd'].unique()))


# In[4]:


print(array_3d[0][0])
print(array_3d.shape)


# In[5]:


bands = 6
pts = 9
dist = 6

for i in range(0, bands):
    for j in range(0, pts):
        ref = array_3d[i][j][0]
        for k in range(0, dist):
            array_3d[i][j][k] = array_3d[i][j][k]-ref

delta = array_3d[0]
theta = array_3d[1]
alpha = array_3d[2]
beta = array_3d[3]
gamma = array_3d[4]
highgamma = array_3d[5]


# In[6]:


ptscols=['aqua','blue','blueviolet','green','orange','red','limegreen','magenta','gold']
ptsids = ['EC133','EC175','EC183','EC186','EC187','EC196','EC219','EC221','EC222']


# In[7]:


fig, axs = plt.subplots(3, 2, figsize=(18, 15))
xsetnw = [-1, 0, 1, 2, 3, 4, 5, 6]
xticks = [' ', 'Referential', '4 mm', '8 mm', '12 mm', '16 mm', '20 mm', ' ']
for i in range(0, 9):  # to 9

    # Delta plot
    axs[0, 0].plot(delta[i], marker='o', color=ptscols[i])
    axs[0, 0].plot(delta[i], color=ptscols[i], label=ptsids[i])
    axs[0, 0].set_xlabel('Bipolar distance')
    axs[0, 0].set_ylabel('Mean log power')
    axs[0, 0].set_title('Delta (2-4Hz)', fontsize=14.5)
    axs[0, 0].set_xlim(-1, 6)
    axs[0, 0].set_ylim(-1.5, 1)  # Normalized y-axis
    axs[0, 0].grid(True)
    axs[0, 0].grid(color='black', linestyle='solid', linewidth=.5)
    axs[0, 0].set_xticks(xsetnw)
    axs[0, 0].set_xticklabels(xticks)

    # Theta plot
    axs[0, 1].plot(theta[i], marker='o', color=ptscols[i])
    axs[0, 1].plot(theta[i], color=ptscols[i], label=ptsids[i])
    axs[0, 1].set_xlabel('Bipolar distance')
    axs[0, 1].set_ylabel('Mean log power')
    axs[0, 1].set_title('Theta (4-8Hz)', fontsize=14.5)
    axs[0, 1].set_xlim(-1, 6)
    axs[0, 1].set_ylim(-1.5, 1)  # Normalized y-axis
    axs[0, 1].grid(True)
    axs[0, 1].grid(color='black', linestyle='solid', linewidth=.5)
    axs[0, 1].set_xticks(xsetnw)
    axs[0, 1].set_xticklabels(xticks)

    # Alpha plot
    axs[1, 0].plot(alpha[i], marker='o', color=ptscols[i])
    axs[1, 0].plot(alpha[i], color=ptscols[i], label=ptsids[i])
    axs[1, 0].set_xlabel('Bipolar distance')
    axs[1, 0].set_ylabel('Mean log power')
    axs[1, 0].set_title('Alpha (8-13Hz)', fontsize=14.5)
    axs[1, 0].set_xlim(-1, 6)
    axs[1, 0].set_ylim(-1.5, 1)  # Normalized y-axis
    axs[1, 0].grid(True)
    axs[1, 0].grid(color='black', linestyle='solid', linewidth=.5)
    axs[1, 0].set_xticks(xsetnw)
    axs[1, 0].set_xticklabels(xticks)

    # Beta plot
    axs[1, 1].plot(beta[i], marker='o', color=ptscols[i])
    axs[1, 1].plot(beta[i], color=ptscols[i], label=ptsids[i])
    axs[1, 1].set_xlabel('Bipolar distance')
    axs[1, 1].set_ylabel('Mean log power')
    axs[1, 1].set_title('Beta (13-25Hz)', fontsize=14.5)
    axs[1, 1].set_xlim(-1, 6)
    axs[1, 1].set_ylim(-1.5, 1)  # Normalized y-axis
    axs[1, 1].grid(True)
    axs[1, 1].grid(color='black', linestyle='solid', linewidth=.5)
    axs[1, 1].set_xticks(xsetnw)
    axs[1, 1].set_xticklabels(xticks)

    # Gamma plot
    axs[2, 0].plot(gamma[i], marker='o', color=ptscols[i])
    axs[2, 0].plot(gamma[i], color=ptscols[i], label=ptsids[i])
    axs[2, 0].set_xlabel('Bipolar distance')
    axs[2, 0].set_ylabel('Mean log power')
    axs[2, 0].set_title('Gamma (25-50Hz)', fontsize=14.5)
    axs[2, 0].set_xlim(-1, 6)
    axs[2, 0].set_ylim(-1.5, 1)  # Normalized y-axis
    axs[2, 0].grid(True)
    axs[2, 0].grid(color='black', linestyle='solid', linewidth=.5)
    axs[2, 0].set_xticks(xsetnw)
    axs[2, 0].set_xticklabels(xticks)

    # High Gamma plot
    axs[2, 1].plot(highgamma[i], marker='o', color=ptscols[i])
    axs[2, 1].plot(highgamma[i], color=ptscols[i], label=ptsids[i])
    axs[2, 1].set_xlabel('Bipolar distance')
    axs[2, 1].set_ylabel('Mean log power')
    axs[2, 1].set_title('High Gamma (50-1200Hz)', fontsize=14.5)
    axs[2, 1].set_xlim(-1, 6)
    axs[2, 1].set_ylim(-1.5, 1)  # Normalized y-axis
    axs[2, 1].grid(True)
    axs[2, 1].grid(color='black', linestyle='solid', linewidth=.5)
    axs[2, 1].set_xticks(xsetnw)
    axs[2, 1].set_xticklabels(xticks)

plt.subplots_adjust(hspace=0.25)
plt.subplots_adjust(wspace=0.15)
#plt.savefig('normalized_power_freqband.png', dpi=300, bbox_inches='tight')
plt.show()


# In[8]:


len(delta), len(delta[0])


# In[9]:


delta


# In[10]:


import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


# In[11]:


def sign_fdr(arr):
    
    ref = []
    
    for i in range(0, len(arr)):
        ref.append(arr[i][0])
        
    comps = np.zeros((5,9))
    for i in range(0, len(arr)):
        for j in range(1, len(arr[i])):
            comps[j-1][i] = arr[i][j]
            
    p_vals = []

    for i in range(0, len(comps)):
        data1 = ref
        data2 = comps[i]
    
        stat, p_value = stats.wilcoxon(data1, data2, alternative='two-sided')
    
        p_vals.append(p_value)
        
    return p_vals
        


# In[12]:


all_ps = []
output_ps = []

all_ps.append(sign_fdr(delta))
all_ps.append(sign_fdr(theta))
all_ps.append(sign_fdr(alpha))
all_ps.append(sign_fdr(beta))
all_ps.append(sign_fdr(gamma))
all_ps.append(sign_fdr(highgamma))

for row in all_ps:
    for c in row:
        output_ps.append(c)


# In[13]:


output_ps


# In[14]:


rejected, p_vals_corrected, _, _ = multipletests(output_ps, alpha=0.05, method='fdr_bh')
print(f'Corrected p-values: {p_vals_corrected}')
print(f'Rejected hypotheses: {rejected}')


# In[15]:


p_vals_corrected.reshape(-1,5)


# In[17]:


p_vals_corrected


# In[18]:


fig, axs = plt.subplots(3, 2, figsize=(18, 15))
xsetnw = [-1, 0, 1, 2, 3, 4, 5, 6]
xticks = [' ', 'Referential', '4 mm', '8 mm', '12 mm', '16 mm', '20 mm', ' ']

freq_bands = [(delta, 'Delta (2-4Hz)'), 
              (theta, 'Theta (4-8Hz)'), 
              (alpha, 'Alpha (8-13Hz)'), 
              (beta, 'Beta (13-25Hz)'), 
              (gamma, 'Gamma (25-50Hz)'), 
              (highgamma, 'High Gamma (50-1200Hz)')]

def plot_band(ax, data, title, i):
    ax.plot(data[i], marker='o', color=ptscols[i])
    ax.plot(data[i], color=ptscols[i], label=ptsids[i])
    ax.set_xlabel('Bipolar distance')
    ax.set_ylabel('Mean log power')
    ax.set_title(title, fontsize=14.5)
    ax.set_xlim(-1, 6)
    ax.set_ylim(-1.5, 1)  
    ax.grid(True)
    ax.grid(color='black', linestyle='solid', linewidth=.5)
    ax.set_xticks(xsetnw)
    ax.set_xticklabels(xticks)

for i in range(0, 9):
    for ax, (data, title) in zip(axs.flat, freq_bands):
        plot_band(ax, data, title, i)

plt.subplots_adjust(hspace=0.25, wspace=0.15)
plt.show()


# In[20]:


import seaborn as sns
import pandas as pd

fig, axs = plt.subplots(2, 3, figsize=(22, 10))
xsetnw = [-1, 0, 1, 2, 3, 4, 5, 6]
xticks = ['4 mm', '8 mm', '12 mm', '16 mm', '20 mm']

freq_bands = [(delta, 'Delta (2-4Hz)'), 
              (theta, 'Theta (4-8Hz)'), 
              (alpha, 'Alpha (8-13Hz)'), 
              (beta, 'Beta (13-25Hz)'), 
              (gamma, 'Gamma (25-50Hz)'), 
              (highgamma, 'High Gamma (50-200Hz)')]

def plot_line_and_violin(ax, data, title, i):

    df = pd.DataFrame(data).melt(var_name='Distance', value_name='Mean log power')
    df = df[df['Distance'] != 0]
    sns.violinplot(x='Distance', y='Mean log power', data=df, ax=ax, color='darkgray', alpha=0.1)

    
    ax.set_xlabel('Bipolar Distance from Referential')
    ax.set_ylabel('Mean log power')
    ax.set_ylim(-1.5, 1)
    ax.set_title(title, fontsize=14.5)
    ax.set_xticks(range(len(xticks)))
    ax.set_xticklabels(xticks)
    ax.grid(True)
    return df

for i in range(0, 9):  # Iterate over subjects
    for ax, (data, title) in zip(axs.flat, freq_bands):
        plot_line_and_violin(ax, data, title, i)

plt.subplots_adjust(hspace=0.28, wspace=0.2)
plt.show()


# In[117]:


# Function to plot both line and violin plots for each subplot
def plot_line_and_violindf(ax, data, title, i):
    
    # Prepare data for violin plot
    df = pd.DataFrame(data).melt(var_name='Distance', value_name='Mean log power')
    df = df[df['Distance'] != 0]
    return df

for i in range(0, 9):  # Iterate over subjects
    for ax, (data, title) in zip(axs.flat, freq_bands):
        df = plot_line_and_violindf(ax, data, title, i)
df

