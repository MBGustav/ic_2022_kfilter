import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv


def read_matrix(fname):
    with open(fname, "r") as f:
        reader  = csv.reader(f, delimiter=",")
        dctread = csv.DictReader(f)
        # HEADER
        row, col = list(map(int, dctread.fieldnames))
        matrix = np.zeros((row, col))
        for i, line in enumerate(reader):
            matrix[i,:] = list(map(float, line))
                        
    return matrix

def plot_graphics(title, mlab_res, c_res):

    dim = min(mlab_res.shape)
    dim1 = min(c_res.shape)
    if(dim != dim1):
        exit

    fig, (ax_n) = plt.subplots(dim1, 1)
    fig.suptitle(title)



    for i in range(dim):
        ax_n[i].plot(mlab_res[:, i], "b")
        ax_n[i].plot(c_res[:, i]   , "r")
        
    plt.show()


euler = read_matrix("DadosZurich/euler.csv")
euler_c = read_matrix("output/euler.csv")
vel_c = read_matrix("output/output_vel.csv")
pos_c = read_matrix("output/output_pos.csv")
qtn_c = read_matrix("output/output_qtn.csv")
rpy_c = read_matrix("output/output_rpy.csv")
bias_gyr_c = read_matrix("output/output_bias_gyr.csv")
bias_acc_c = read_matrix("output/output_bias_acc.csv")

vel_m = read_matrix("output/mlab_vel.csv")
pos_m = read_matrix("output/mlab_pos.csv")
qtn_m = read_matrix("output/mlab_qtn.csv")
rpy_m = read_matrix("output/mlab_rpy.csv")
bias_gyr_m = read_matrix("output/mlab_bias_gyr.csv")
bias_acc_m = read_matrix("output/mlab_bias_acc.csv")


plot_graphics("Roll, pitch, yaw", euler,euler_c)
plot_graphics("Velocity(x,y,z)", vel_m,vel_c)
plot_graphics("Position(x,y,z)", pos_m,pos_c)
plot_graphics("Quarternions", qtn_m,qtn_c)
# plot_graphics("RPY", rpy_m,rpy_c)
plot_graphics("bias acc", bias_acc_m,bias_acc_c)
plot_graphics("bias gyro", bias_gyr_m,bias_gyr_c)