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
pmm = load_trajectory_samples_pmm(trajectory_file_pmm)


fig, axs = plt.subplots(8)
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
    

min_distances = []
min_distances_pos = []

    
md_max = 0;
for md in min_distances:
    if md > md_max:
        md_max = md

md_pos_max = 0;
for mdp in min_distances_pos:
    if mdp > md_pos_max:
        md_pos_max = mdp

print("md_max ",md_max)
print("md_pos_max ",md_pos_max)

plt.legend()

plt.show()  
