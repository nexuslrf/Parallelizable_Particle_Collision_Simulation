#!/usr/bin/env python3

# Make sure to install all the necessary python libraries via pip
# and ImageMagick via your OS package manager

# Authors:
# LEE YONG JIE, RICHARD
# KEVEN LOO YUQUAN


import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import pi
from time import time
import sys

if len(sys.argv) > 1:
    filepath = sys.argv[1]
else:
    filepath = input(
        'Enter the input/output file you would like to visualize (xxx.in & xxx.out): '
    )
# filepath = 'testAnim'

_SLOWDOWN_FACTOR = 1

startTime = time()

for extension in ('.in', '.out'):
    try:
        open(filepath + extension).close()
    except FileNotFoundError:
        print('{}{} not found, exiting...'.format(filepath, extension))
        exit(0)

dataDict = {}

with open(filepath + '.in', 'r') as particleIn:
    N = int(particleIn.readline().strip())
    L = float(particleIn.readline().strip())
    r = float(particleIn.readline().strip())
    steps = int(particleIn.readline().strip())

steps *= _SLOWDOWN_FACTOR
maxSpeed = L / 4 / _SLOWDOWN_FACTOR

with open(filepath + '.out', 'r') as particleOut:
    for line in particleOut.readlines():
        data = line.strip().split()

        step, index, x, y, v_x, v_y, *others = data

        step = int(step)
        index = int(index)
        x = float(x)
        y = float(y)
        v_x = float(v_x)
        v_y = float(v_y)
        v = (v_x**2 + v_y**2)**0.5

        if step not in dataDict:
            dataDict[step] = {}
        dataDict[step][index] = (x, y, v)

for i in range(steps + 1):
    if i not in dataDict:
        print('WARNING: Timestep {} not in output'.format(i))
    if i in dataDict and len(dataDict[i]) != N:
        print('WARNING: Unequal number of particles (detected at timestep {})'.
              format(i))

fig, ax = plt.subplots(figsize=(8, 8))
time_text = ax.text(0.8,
                    1.05,
                    '',
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    size=14)


def init():
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    return []


circlePatch = []


def plotPoints(step):
    circles = []
    time_text.set_text('Frame: {}'.format(step))
    circles.append(time_text)
    for patch in circlePatch:
        patch.remove()
    circlePatch.clear()
    for v in dataDict[step].values():
        scale = max(min(v[2] / maxSpeed, 1), 0)
        newPatch = plt.Circle(v[:2], r, color=(scale, 0, 1 - scale))
        circlePatch.append(newPatch)
        circles.append(ax.add_patch(newPatch))
    # scatter.set_offsets(list(dataDict[step].values()))
    return circles


ani = FuncAnimation(fig,
                    plotPoints,
                    frames=list(range(0, steps + 1)),
                    init_func=init,
                    blit=True,
                    interval=10)
ani.save('{}.gif'.format(filepath),
         writer='imagemagick',
         fps=2 * _SLOWDOWN_FACTOR)

endTime = time()

print('Written to {}.gif!\n----------'.format(filepath))
print('Time taken: {:.3f}s'.format(endTime - startTime))
