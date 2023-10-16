import matplotlib.pyplot as plt
import numpy as np

def get_station_colors():
    cmx = plt.get_cmap('inferno')
    cols = [cmx(i) for i in np.linspace(0.1, .85, 5)]
    #colours = ['darkorchid', '#377eb8', '#4dac26', '#ca0020'] # [st, wb, o8, tr]
    #cols = ['#4dac26', '#377eb8', '#e66101', 'darkorchid', 'black']
    cols = ['#4dac26', '#377eb8', 'orange', 'darkorchid', 'black']
    return {'on': cols[0], 'wb': cols[1], 'tr': cols[2], 'st': cols[3], 'fast': cols[4]}
