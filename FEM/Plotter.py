import matplotlib.pyplot as plt

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
