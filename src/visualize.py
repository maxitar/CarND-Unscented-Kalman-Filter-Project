import numpy as np
import matplotlib.pyplot as plt
import csv
import os

sensor_type = []
times = []
data = []
with open('output.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter = '\t')
    names = reader.__next__()
    for line in reader:
        sensor_type.append(line[0])
        times.append(line[1])
        data.append(line[2:])

names.pop(0)
names.pop(0)
sensor_type = np.array(sensor_type)
times = np.array(times, dtype=float)
data = np.array(data, dtype=float)

titles = ['px','py','vx','vy','yaw', 'yaw_rate']
units = ['m','m','m/s','m/s','rad','rad/s']
fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(20,15))
current_datum = 0
for row in range(3):
    for col in range(2):
        diff = data[:,current_datum] - data[:,current_datum+1]
        rmse = np.sqrt(np.mean(diff*diff))
        print("{} RMSE: {:.3f}".format(titles[row*2+col],rmse))
        ax[row, col].set_title(titles[row*2+col])
        ax[row, col].set_xlabel('s')
        ax[row, col].set_ylabel(units[row*2+col])
        ax[row, col].plot(times, data[:,current_datum], '-.', label = names[current_datum])
        ax[row, col].plot(times, data[:,current_datum+1], label = names[current_datum + 1])
        ax[row, col].legend()
        current_datum += 2

fig.savefig('plots.png', bbox_inches='tight')

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(20,15))
lidar_nis = data[sensor_type=='L', 12]
num_sample = lidar_nis.shape[0]
x_axis = np.arange(num_sample)
ax[0].set_title("LIDAR NIS")
ax[0].plot(x_axis, lidar_nis)
ax[0].plot([0, num_sample], [5.991, 5.991])

radar_nis = data[sensor_type=='R', 12]
radar_nis = radar_nis[1:]
num_sample = radar_nis.shape[0]
x_axis = np.arange(num_sample)
ax[1].set_title("RADAR NIS")
ax[1].plot(x_axis, radar_nis)
ax[1].plot([0, num_sample], [7.815, 7.815])

fig.savefig('nis.png', bbox_inches='tight')
