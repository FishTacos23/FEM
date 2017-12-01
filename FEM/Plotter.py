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


def animated_plot(x, dx, y, dy, z, dz):

    # Attaching 3D axis to the figure
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    line, = ax.plot(z, y, zs=x)

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, 0.0, '', transform=ax.transAxes)

    # Setting the axes properties
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_title('3D Test')

    intervals = 50

    dt = 0.05
    t = np.arange(0.0, 20, dt)

    def animate(i):
        factor = math.sin(2. * math.pi * float(i) / float(intervals))
        this_x = []
        this_y = []
        this_z = []

        for j in xrange(len(x)):
            this_x.append(x[j] + dx[j] * factor)
            this_y.append(y[j] + dy[j] * factor)
            this_z.append(z[j] + dz[j] * factor)

        line.set_data(this_z, this_y)
        time_text.set_text(time_template % (i*dt))
        return line, time_text

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    # Creating the Animation object
    ani = animation.FuncAnimation(fig, animate, frames=100, interval=50, blit=True, init_func=init)
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

    for i, p in enumerate(mx):
        plt.plot(n, p, color=colors[i], marker='+', label=labels[i])

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
