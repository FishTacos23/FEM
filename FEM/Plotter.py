import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import math

colors = ['r', 'g', 'b', 'orange', 'black', 'cyan', 'magenta', 'yellow', 'darkkhaki', 'darksalmon']


def plot_graphs(x, y, title, p=None):

    for i in xrange(len(x)):
        for j in xrange(len(x[i])):
            plt.plot(x[i][j], y[i][j], color=colors[i])

    if p is not None:
        plt.plot(p[1], p[0], 'bo')

    plt.xlabel('x')
    plt.ylabel('u')
    plt.title(title)
    plt.show()


def plot_deflections(n, mx, exact, title, labels):

    for i, p in enumerate(mx):
        plt.plot(n, p, color=colors[i], marker='+', label=labels[i])

    if exact is not None:
        plt.plot((0, max(n)), (exact, exact), color='b', label='Exact')

    plt.xlabel('Number of Nodes')
    plt.ylabel('Tip Deflection')
    plt.title(title)
    plt.show()


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
