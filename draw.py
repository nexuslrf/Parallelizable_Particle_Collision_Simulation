'''
if inputs is in print mode,
you can use this scrip to visualize collision
'''

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm
from matplotlib.patches import Ellipse, Circle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', type=str, default='inputs.txt')
parser.add_argument('--output_file', type=str, default='outputs.txt')
parser.add_argument('--fig_name', type=str, default='colli_demo.png')
parser.add_argument('--wide', type=int, default=3)
parser.add_argument('--height', type=int, default=3)
args = parser.parse_args()

params = [eval(i) if i.isdigit() else i for i in open(args.input_file).readline().split()]
num_p = params[0]
l = params[1]
r = params[2]
steps = params[3]
# print(params)
h = args.wide
w = args.height
particles = np.array([
    [eval(i) for i in line.split()[1:]]
    for line in open(args.output_file).readlines()[:num_p*h*w]
])

norm = mpl.colors.Normalize(vmin=0, vmax=num_p)
cmap = cm.rainbow
mps = cm.ScalarMappable(norm=norm, cmap=cmap)
fig, ax = plt.subplots(h,w, figsize=(3*h,3*w))
for i in range(h):
    for j in range(w):
        ax[i][j].set_xlim(0,l)
        ax[i][j].set_ylim(0,l)
        sca = ax[i][j].scatter(particles[(h*i+j)*num_p:(h*i+j+1)*num_p,1], 
                         particles[(h*i+j)*num_p:(h*i+j+1)*num_p,2], 
                         s=2, c=particles[(h*i+j)*num_p:(h*i+j+1)*num_p,0], 
                         cmap=cmap)
        for k in range((h*i+j)*num_p, (h*i+j+1)*num_p):
            cir = Circle(xy=(particles[k][1], particles[k][2]), alpha=0.4, 
                radius=r, color=mps.to_rgba(particles[k][0]))
            ax[i][j].add_patch(cir)
            ax[i][j].quiver(particles[k,1],particles[k,2], 
                        particles[k,3], particles[k,4], 
                        color=mps.to_rgba(particles[k,0]))
        ax[i][j].set_title("Step {}".format(h*i+j))

fig.tight_layout()
fig.subplots_adjust(right=0.8,bottom=0.2)
cbar_ax2 = fig.add_axes([0.85, 0.2, 0.05, 0.75])
cbar = fig.colorbar(sca, cax=cbar_ax2, ticks=[0,4,9,14,19])
cbar.ax.set_ylabel('Particle Index', rotation=270)
# plt.tight_layout()
# plt.show()
plt.savefig(args.fig_name)