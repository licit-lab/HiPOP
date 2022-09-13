from hipop.cpp.graph import OrientedGraph

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def render_oriented_graph(ax, G: OrientedGraph, color='black', linkwidth=1, nodesize=2, node_label=True, show_length=False, cmap=plt.cm.jet):
    x, y = zip(*[n.position for n in G.nodes.values()])
    ax.plot(x, y, 'o', markerfacecolor='white', markeredgecolor=color, fillstyle='full', markersize=nodesize)

    lines = []
    for lid, link in G.links.items():
        lines.append([G.nodes[link.upstream].position, G.nodes[link.downstream].position])
    line_segment = LineCollection(lines, linestyles='solid', linewidths=linkwidth, cmap=cmap)
    ax.add_collection(line_segment)





if __name__ == "__main__":
    from hipop.graph import generate_manhattan

    G = generate_manhattan(3, 10)

    fig, ax = plt.subplots()
    render_oriented_graph(ax, G)
    plt.show()