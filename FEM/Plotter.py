import matplotlib.pyplot as plt

colors = ['r', 'g', 'b', 'orange']


def plot_graphs(x, y):

    for i in xrange(len(x)):
        for j in xrange(len(x[i])):
            plt.plot(x[i][j], y[i][j], color=colors[i])

    plt.show()
