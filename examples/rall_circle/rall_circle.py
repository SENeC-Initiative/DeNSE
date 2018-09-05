import dense as ds
import sys
import matplotlib.pyplot as plt
import matplotlib
import copy

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def rall_circle(pop,gid):
    ax = plt.subplot(221)
    ax.set_title("Dendritic density field")
    bx = plt.subplot(222)
    bx.set_title("Axon density field")
    cx = plt.subplot(223, projection='polar')
    cx.set_title("Sholl Analysis: dendrites")
    # label_position=cx.get_rlabel_position()
    label = r'soma distance\mu m'
    import numpy as np
    dx = plt.subplot(224, projection='polar')
    dx.set_title("Sholl Analysis: axon")
    # label_position=cx.get_rlabel_position()
    for ng_neuron in pop[10:15]:
        for dendrite in ng_neuron.dendrites:
            # ax.scatter(dendrite.xy[:,0] - ng_neuron.position[0],
                        # dendrite.xy[:,1] - ng_neuron.position[1], c='b',alpha=0.003,
                       # edgecolors=None)
            cx.scatter(dendrite.branching_points[:,0] - ng_neuron.position[0],
                        dendrite.branching_points[:,1] - ng_neuron.position[1], edgecolors=None)

        # bx.scatter(axon.xy[:,0] - ng_neuron.position[0],
                     # axon.xy[:,1] - ng_neuron.position[1]
                     # ,' c='r', alpha=0.0051,edgecolors=None)
        axon = ng_neuron.axon
        dx.scatter(axon.branching_points[:,0] - ng_neuron.position[0],
                    axon.branching_points[:,1] - ng_neuron.position[1]
                     ,edgecolors=None)
    from matplotlib.colors import LogNorm
    my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
    my_cmap.set_bad((0,0,0))

    _axon_ = pop.axon_all_points(center_zero=True)
    _dendrites_ = pop.dendrites_all_points(center_zero=True)
    ax.hist2d(_dendrites_[:,0], _dendrites_[:,1],bins=100,
            norm=LogNorm(),range=[[-5000,5000],[-5000,5000]], cmap=my_cmap)
    bx.hist2d(_axon_[:,0],_axon_[:,1], bins=100, norm=LogNorm(),\
              range=[[-5000,5000],[-5000,5000]], cmap=my_cmap)
    ax.set_xlabel("X")
    bx.set_xlabel("X")
    ax.set_ylabel("Y")
    bx.set_ylabel("Y")

    cx.text(np.radians(-30),cx.get_rmax(),label,
        rotation=0.,ha='left',va='top')
    dx.text(np.radians(-30),dx.get_rmax(),label,
        rotation=0.,ha='left',va='top')
    plt.tight_layout()
    plt.show()
    # plt.savefig("sholl_analysis.png", format='png',ppi=300)
        # axes.scatter(0,0,c='k')

if __name__ =="__main__":
    pop = ds.Population(sys.argv[1])
    import argparse
    parser = argparse.ArgumentParser(description='Sholl analysis')
    parser.add_argument('--culture',  type=str,
                        help='folder with swc files')
    parser.add_argument('--no_import', action='store_false', default=True,
                        help='do not import the population from file')
    args = parser.parse_args()

    if args.no_import:
        for n in range(6):
            swc_culture =ds.NeuronsFromSimulation(args.culture+str(n))
            pop.info = swc_culture["info"]
            pop.add_swc_population(swc_culture['neurons'])


    print("population imported: with {} neurons".format(len(pop.gids)))
    gids = list(pop.gids)
    rall_circle(pop, gids)


