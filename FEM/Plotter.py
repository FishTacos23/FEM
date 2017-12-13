import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import math

colors = ['r', 'g', 'b', 'orange', 'black', 'cyan', 'magenta', 'yellow', 'darkkhaki', 'darksalmon']


def plot_modes(x, y_list):
    for y in xrange(len(y_list)):
        for e in xrange(len(x)):
            plt.plot(x[e], y_list[y][e], color=colors[y])
    plt.show()


def plot_errors(x, y):
    for i in xrange(len(x)):
        plt.plot(x[i], y[i])
    plt.show()


def animation_plot(x, y):

    dt = 0.05

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(0, 1), ylim=(-.5, .5))
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def animate(i):
        this_y = y*math.cos(2.*math.pi*i/100)

        line.set_data(x, this_y)
        time_text.set_text(time_template % (i*dt))
        return line, time_text

    ani = animation.FuncAnimation(fig, animate, np.arange(1, 100), interval=50, blit=True, init_func=init)
    # ani.save('double_pendulum.mp4', fps=15)
    plt.show()
