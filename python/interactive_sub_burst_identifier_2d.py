import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


def ds(array, factor=2, axis=0):
    """Downsample using the mean along a given axis"""
    if factor > 1:
        if axis < 0:
            axis += array.ndim

        # Trim the array if necessary.
        axis_len = array.shape[axis]
        if axis_len % factor != 0:
            print("The array axis with length %s is trimmed to fit the downsample factor "%(axis_len))
            slc = [slice(None)] * array.ndim
            slc[axis] = slice(-(axis_len % factor))
            array = array[tuple(slc)]

        new_shape = list(array.shape)
        new_shape[axis] //= factor
        new_shape.insert(axis+1, factor)
        array = array.reshape(new_shape).mean(axis=axis+1)
    return array


class identify_bursts(object):
    def __init__(self, arr, tsamp, title='Interactive plotter', block=True):
        self.peak_times = []
        self.peak_freqs = []
        self.peak_amps = []
        self.tsamp = tsamp
        self.tavg = 1
        self.favg = 1
        self.arr = arr
        self.vmin_fac = 1.

        profile = np.sum(arr, axis=0)
        spectrum = np.sum(arr, axis=1)

        # Plot the whole thing
        rows = 2
        cols = 2
        gs = gridspec.GridSpec(rows, cols, wspace=0., hspace=0.,
                               height_ratios=[0.5,]*(rows-1) + [2,],
                               width_ratios=[5,] + [1,]*(cols-1))

        figs = plt.figure(title, figsize=(8, 5))
        plt.title(title)
        self.ax_data = plt.subplot(gs[2])  # dynamic spectrum
        self.ax_ts = plt.subplot(gs[0], sharex=self.ax_data)  # time series
        self.ax_spec = plt.subplot(gs[-1], sharey=self.ax_data)  # spectrum
        self.canvas = self.ax_data.figure.canvas

        self.ts_plot, = self.ax_ts.plot(profile, 'k-', alpha=1.0)
        y_range = profile.max() - profile.min()
        self.ax_ts.set_ylim(profile.min() - y_range * 0.1, profile.max() + y_range * 0.1)

        self.spec_plot, = self.ax_spec.step(spectrum, np.arange(spectrum.shape[0]), 'k-')
        y_range = spectrum.max() - spectrum.min()
        self.ax_spec.set_xlim(spectrum.min() - y_range * 0.1, spectrum.max() + y_range * 0.1)

        self.data_plot = self.ax_data.imshow(arr, vmin=-1*arr.std(), vmax=7*arr.std(),
                                             origin='lower', aspect='auto',interpolation='nearest')

        plt.colorbar(self.data_plot)

        plt.tight_layout()

        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.keyPress = self.canvas.mpl_connect('key_press_event', self.onKeyPress)

        plt.show(block=block)

    def onpress(self, event):
        tb = plt.get_current_fig_manager().toolbar
        if tb.mode == '':
            y1 = int(round(event.ydata))
            x1 = event.xdata
            self.peak_times.append(x1)
            self.peak_freqs.append(y1)
            meanval = np.mean(self.arr[y1-5:y1+6, int(x1)-5:int(x1)+6])
            self.peak_amps.append(meanval)  # correct for downsampling
            self.ax_data.scatter(x1, y1, lw=1, color='r', marker='x', s=100,
                                 zorder=10)
            print(x1, y1, meanval)
            plt.draw()
    def onKeyPress(self, event):
        if event.key == 't':
            self.tavg *= 2
        elif event.key == 'F':
            self.favg *= 2
        elif event.key == 'u':
            if self.tavg > 1:
                self.tavg //= 2
        elif event.key == 'j':
            if self.favg > 1:
                self.favg //= 2
        elif event.key == 'e':
            self.vmin_fac -= .1
            if self.vmin_fac < .1:
                self.vmin_fac = 1.

        arr = ds(ds(self.arr, factor=self.favg), factor=self.tavg, axis=1)
        if event.key in 'tFuj':
            profile = np.sum(arr, axis=0)
            spectrum = np.sum(arr, axis=1)

            # Replot.
            self.ts_plot.set_data(np.arange(profile.shape[0]*self.tavg, step=self.tavg), profile)
            y_range = profile.max() - profile.min()
            self.ax_ts.set_ylim(profile.min() - y_range * 0.1, profile.max() + y_range * 0.1)

            self.spec_plot.set_data(spectrum, np.arange(spectrum.shape[0]*self.favg, step=self.favg))
            y_range = spectrum.max() - spectrum.min()
            self.ax_spec.set_xlim(spectrum.min() - y_range * 0.1, spectrum.max() + y_range * 0.1)

            self.data_plot.set_data(arr)

        self.data_plot.set_clim(vmin=self.vmin_fac*arr.min(), vmax=arr.max())
        plt.draw()
