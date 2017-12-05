import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import math
import numpy as np

colors = ['r', 'g', 'b', 'orange', 'black']


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


def new_plot(x, dx, y, dy, z, dz):
    f = 100
    tx = []
    tz = []

    for ii in xrange(len(x)):
        tx.append(x[ii]+dx[ii]*f)
        tz.append(1.)

    plt.plot(x, z, 'r+')
    plt.plot(tx, tz, 'bv')
    plt.show()


def animated_plot(x, dx, y, dy):

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
    ax.grid()

    line, = ax.plot([], [], 'r+', lw=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    ax.set_xlim([(min(x)-min(dx))*1.1, (max(x)+max(dx))*1.1])

    # Setting the axes properties
    intervals = 50

    dt = 0.05
    t = np.arange(0.0, 20, dt)

    def animate(q):
        factor = math.sin(2. * math.pi * float(q) / float(intervals))
        this_x = []
        this_y = []

        for j in xrange(len(x)):
            this_x.append(x[j] + dx[j] * factor)
            this_y.append(y[j] + dy[j] * factor)

        line.set_data(this_x, this_y)
        time_text.set_text(time_template % (q*dt))
        return line, time_text

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    # Creating the Animation object
    ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)), interval=intervals, blit=True, init_func=init)
    plt.show()


def gen_lines(dx, dy, dz, x, y, z, intervals):
    data = []
    for i in xrange(intervals):
        line = [[], [], []]
        factor = math.sin(2.*math.pi*float(i)/float(intervals))
        for j in xrange(len(x)):
            line[0].append(x+dx*factor)
            line[1].append(y+dy*factor)
            line[2].append(z+dz*factor)
        data.append(line)
    return data


def update_animation(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines


def plot_deflections(n, mx, exact, title, labels):

    print exact

    for i, p in enumerate(mx):
        plt.plot(n[i], p, color=colors[i], marker='o', label=labels[i])

    if exact is not None:
        plt.plot((0, max(n)), (exact, exact), color='b', label='Exact')

    plt.xlabel('Number of Nodes')
    plt.ylabel('Tip Deflection')
    plt.title(title)
    plt.show()


def plt_error(h, e, title):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for i in xrange(len(e)):
        line, = ax.plot(h, e[i], color=colors[i], lw=2)

    ax.set_yscale('log')
    ax.set_xscale('log')
    fig.suptitle(title)

    plt.show()
