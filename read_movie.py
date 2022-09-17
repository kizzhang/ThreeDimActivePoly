# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 16:59:45 2022

@author: Zhiyu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.animation as animation

time = 10
df = pd.read_csv('S:\Programs\Active Polymer\ABP polymer\walk_trajectory.txt', header = None,index_col= False)
num = df.iloc[0].to_numpy(dtype ='int')[0]
x = np.zeros((time,num,3))
x_ani = []

def update_plot(num, walks, lines):
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[num, :2].T)
        line.set_3d_properties(walk[num, 2])
    return lines

for i in range(time):
    for j in range(num):
        x[i][j][0] = float(df.iloc[i+i*num+j+3].str.split(pat = '\t',).to_numpy()[0][1])
        x[i][j][1] = float(df.iloc[i+i*num+j+3].str.split(pat = '\t',).to_numpy()[0][2])
        x[i][j][2] = float(df.iloc[i+i*num+j+3].str.split(pat = '\t',).to_numpy()[0][3])
    x_ani.append(x[i])
    #print(x_ani)


## MAKE ANIMATION 

#print(x)
# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# Create lines initially without data
lines = [ax.plot([], [], [],'o')[0] for _ in x_ani]

# Setting the axes properties
ax.set(xlim3d=(min(x[:,:,0].flatten()), max(x[:,:,0].flatten())), xlabel='X')
ax.set(ylim3d=(min(x[:,:,1].flatten()), max(x[:,:,1].flatten())), ylabel='Y')
ax.set(zlim3d=(min(x[:,:,2].flatten()), max(x[:,:,2].flatten())), zlabel='Z')

# Creating the Animation object
ani = animation.FuncAnimation(fig, update_plot, time, fargs=(x_ani, lines), interval=300)

plt.show()

# writer = animation.writers['ffmpeg'](fps=10)

# ani.save('demo.mp4',writer=writer,dpi=100)
