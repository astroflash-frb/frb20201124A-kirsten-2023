"""
To plot the time series of a dataset to identify the sub-burst peaks
"""
import matplotlib.pyplot as plt


class identify_bursts(object):
    def __init__(self, arr, title='Identify bursts', pltrange=None):
        self.peak_times = []
        self.peak_amps = []
        self.ranges = []
        self.lines = []

        profile = arr

        fig = plt.figure(title, figsize=(8, 5))
        plt.title(title)

        self.ax_ts = plt.subplot(111)
        self.canvas = self.ax_ts.figure.canvas

        self.ts_plot, = self.ax_ts.plot(profile, 'k-', alpha=1.0)
        y_range = profile.max() - profile.min()
        self.ax_ts.set_ylim(profile.min() - y_range * 0.1, profile.max() + y_range * 0.1)

        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.keyPress = self.canvas.mpl_connect('key_press_event', self.onKeyPress)
        self.keyRelease = self.canvas.mpl_connect('key_release_event', self.onKeyRelease)
        self.i = False
        if pltrange is not None:
            plt.xlim(pltrange)
        plt.show()

    def onpress(self, event):
        if self.i:
            tb = plt.get_current_fig_manager().toolbar
            if tb.mode == '':
                y1 = event.ydata
                x1 = event.xdata
                self.peak_times.append(x1)
                self.peak_amps.append(y1)

                lx = self.ax_ts.axvline(x1, lw=1, color='yellow', zorder=10)
                ly = self.ax_ts.axhline(y1, lw=1, color='yellow', zorder=10)
                #self.ax_ts.scatter(x1, y1, lw=1, color='r', marker='x', s=100,zorder=10)
                print(x1, y1)
                plt.draw()
        else:
            tb = plt.get_current_fig_manager().toolbar
            if tb.mode == '':
                x1 = event.xdata
                self.ranges.append(int(x1))
                l = self.ax_ts.axvline(x1, lw=1, color='purple',zorder=10)
                self.lines.append(l)
                plt.draw()

    def onKeyPress(self, event):
        if event.key == 'i':
            self.i = True
        if event.key == 'u':
            del self.ranges[-1]
            print(f'Removed last entry from list')
            # TODO: figure out how to remove the corresponding line and replot...

    def onKeyRelease(self, event):
        if event.key == 'i':
            self.i = False
