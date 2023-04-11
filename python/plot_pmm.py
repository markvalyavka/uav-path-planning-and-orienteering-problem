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
    # print("load_trajectory_samples ",file)
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
    return np.array(states)


def get_trajectory_positions(config_file):
    with open(config_file) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
        locations = config['locations']
        end = config['end']['position']
        start = config['start']['position']
        rewards = config['rewards']
        return rewards, np.array([start, *locations, end])


def plot_3d_positions_graph(pmm, pmm_equidistant):
    fig2 = plt.figure(figsize=(12, 7))
    ax3d = fig2.add_subplot(111, projection='3d', computed_zorder=False)
    fig2.text(0, 1, f"Reward: 350", transform=ax3d.transAxes, fontsize=10,
              verticalalignment='center', bbox=dict(facecolor='#f7f7f7', edgecolor='#333333', linewidth=1.5, boxstyle='round,pad=0.5'))

    velocity_norms = np.sqrt(pmm[:, 4] * pmm[:, 4] + pmm[:, 5] * pmm[:, 5] + pmm[:, 6] * pmm[:, 6])
    # for i, v in enumerate(velocity_norms):
    #     print(i, v)
    # ax3d.plot(pmm[:, 1], pmm[:, 2], pmm[:, 3])
    velocities_plot = ax3d.scatter(pmm[:, 1], pmm[:, 2], pmm[:, 3], c=velocity_norms, cmap='jet', s=0.35, zorder=-1)

    # ax3d.plot(pmm_equidistant[:, 1], pmm_equidistant[:, 2], pmm_equidistant[:, 3], '.k')

    rewards, locations = get_trajectory_positions('../new_config.yaml')
    # print(locations)
    p = ax3d.scatter(locations[1:len(locations)-1, 0], locations[1:len(locations)-1, 1], locations[1:len(locations)-1, 2], c=rewards, cmap='viridis_r', s=[30]*(len(locations)-2), alpha=1)
    ax3d.scatter(locations[0, 0], locations[0, 1], locations[0, 2], c='r', s=[30], alpha=1)
    ax3d.scatter(locations[len(locations)-1, 0], locations[len(locations)-1, 1], locations[len(locations)-1, 2], c='r', s=[30], alpha=1)

    ax3d.set_xlabel('x in m', labelpad=10)
    ax3d.set_ylabel('y in m', labelpad=10)
    ax3d.set_zlabel('z in m', labelpad=10)
    fig2.colorbar(p, label="rewards", ticks=[int(max(rewards)/4), int(max(rewards)/2), int(max(rewards)*3/4), int(max(rewards))], shrink=0.7)
    fig2.colorbar(velocities_plot, label="velocity m/s", pad=0.15,
                  ticks=[int(max(velocity_norms) / 4), int(max(velocity_norms) / 2), int(max(velocity_norms) * 3 / 4), int(max(velocity_norms))], shrink=0.7)
    fig2.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.05, hspace=0.2)
    plt.legend()
    plt.show()


# fig, axs = plt.subplots(7)
# ax = fig.add_subplot(111, projection='3d')

# colors = ['k','c','r','b']
# max_ax=0
# min_ax=0

#print("sequence",sequence)
#lc = Line3DCollection(sequence, colors = 'red')
#ax.add_collection(lc)

# for i in range(6):
#     label_pmm = 'pmm p(%i)'%(i)
#     axs[i].plot(pmm[:,0]/pmm[-1,0],pmm[:,i+1],'g',label=label_pmm)
#     if i < 3:
#         axs[i].plot(pmm_equidistant[:,0]/pmm_equidistant[-1,0],pmm_equidistant[:,i+1],'.k',label=label_pmm)
#
# for i in range(1,pmm_equidistant.shape[0]):
#     d = np.linalg.norm(pmm_equidistant[i,1:4]-pmm_equidistant[i-1,1:4])
#     axs[6].plot(pmm_equidistant[i,0]/pmm_equidistant[-1,0],d,'.')


if __name__ == "__main__":
    trajectory_pmm_file = '../samples_pmm.csv'
    trajectory_pmm_equidistant_file = '../samples_equidistant.csv'
    pmm = load_trajectory_samples_pmm(trajectory_pmm_file)
    # pmm_equidistant = load_trajectory_samples_pmm(trajectory_pmm_equidistant_file)

    plot_3d_positions_graph(pmm, None)

