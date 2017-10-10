import matplotlib.pyplot as plt

colors = ['r', 'g', 'b', 'orange']


def plot_graphs(x, y, title):

    for i in xrange(len(x)):
        for j in xrange(len(x[i])):
            plt.plot(x[i][j], y[i][j], color=colors[i])

    plt.xlabel('x')
    plt.ylabel('u')
    plt.title(title)
    plt.show()


def plt_error(h, e, labels, title):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for i in xrange(len(e)):
        line, = ax.plot(h, e[i], color=colors[i], lw=2)

    ax.set_yscale('log')
    ax.set_xscale('log')
    fig.suptitle(title)

    plt.show()
