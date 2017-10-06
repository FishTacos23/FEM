import matplotlib.pyplot as plt


def plot_graphs(x, y):

    for i in xrange(len(x)):
        plt.plot(x[i], y[i])

    plt.show()
