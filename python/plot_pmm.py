#!/usr/bin/env python3

import random, sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import plot
from fileinput import filename
from cmath import acos, pi, cos, sin
from math import atan2, sqrt
import csv, copy
from mpl_toolkits.mplot3d import axes3d
import yaml


def load_trajectory_samples_pmm(file):
    print("load_trajectory_samples ",file)
    states = []
    with open(file, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)
        last_pos = None
        for row in csvreader:
            col = []
            for c in row:
                col.append(float(c))
            state = [col[0],col[1],col[2],col[3],col[8],col[9],col[10]]
            states.append(state)    
            
    #print(edges)
    return np.array(states)


gates = None
end = None
start = None
with open('../config.yaml') as file: 
    config = yaml.load(file,Loader=yaml.FullLoader)
    print('start',config['start']['position']) 
    print('end',config['end']['position'])
    print("gates",config['gates'])
    gates = config['gates']
    end = config['end']['position']
    start = config['start']['position']


trajectory_file_pmm = '../samples_pmm.csv'
trajectory_file_pmm_equidistant = '../samples_equidistant.csv'
pmm = load_trajectory_samples_pmm(trajectory_file_pmm)
pmm_equidistant = load_trajectory_samples_pmm(trajectory_file_pmm_equidistant)


fig, axs = plt.subplots(7)
# ax = fig.add_subplot(111, projection='3d')

colors = ['k','c','r','b']
max_ax=0
min_ax=0
        
#print("sequence",sequence)
#lc = Line3DCollection(sequence, colors = 'red')
#ax.add_collection(lc)

for i in range(6):
    label_pmm = 'pmm p(%i)'%(i)
    axs[i].plot(pmm[:,0]/pmm[-1,0],pmm[:,i+1],'g',label=label_pmm)
    if i < 3:
        axs[i].plot(pmm_equidistant[:,0]/pmm_equidistant[-1,0],pmm_equidistant[:,i+1],'.k',label=label_pmm)
    
for i in range(1,pmm_equidistant.shape[0]):
    d = np.linalg.norm(pmm_equidistant[i,1:4]-pmm_equidistant[i-1,1:4])
    axs[6].plot(pmm_equidistant[i,0]/pmm_equidistant[-1,0],d,'.')

fig2 = plt.figure()
ax3d = fig2.add_subplot(111, projection='3d')
ax3d.plot(pmm[:,1],pmm[:,2],pmm[:,3])
ax3d.plot(pmm_equidistant[:,1],pmm_equidistant[:,2],pmm_equidistant[:,3],'.k')

plt.legend()

plt.show()  
