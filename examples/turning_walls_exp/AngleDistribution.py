import NetGrowth
import warnings
from NetGrowthRwBenchmark import AnalyseNetgrowthRW, CleanFolder
import numpy as np
import matplotlib.pyplot as plt
import os

n_samples = 200
angle = 90 # 20, 40, 90
culture_file = "../culture/angle"+str(angle)+".svg"
# this is linked to the file .svg
end_of_the_tube = 300
sim_length = 320
cavity_x_max = 600
neuron_x_max = 50
filopodia_length = 10. ##PUT the DOT!
save_path ="renaud_angle"+str(angle)+"_fil"+str(filopodia_length)


def RunNetGrowth(n_samples, sim_length, n_procs, neuron_params,
                 save_path="tmp_measure", plot=False):
    """
    Run NetGrowth simulation
    """
    kernel = {"seeds": [33, 57, 19, 37, 79, 87, 11][:n_procs],
              "num_local_threads": n_procs,
              "resolution": 1.}
    experiment_params = {}
    experiment_params["num_neurons"] = n_samples
    np.random.seed(kernel['seeds'])
    NetGrowth.SetKernelStatus(kernel, NetGrowth.GenerateSimulationID())
    culture = NetGrowth.CreateEnvironment(
        culture_file, min_x=0, max_x=cavity_x_max)
    pos_left = culture.seed_neurons(
        neurons=experiment_params["num_neurons"], xmax=neuron_x_max,
        soma_radius=1.)

    neuron_params['growth_cone_model'] = 'random_walk'
    neuron_params['position'] = pos_left

    gids = None
    gids = NetGrowth.CreateNeurons(experiment_params["num_neurons"],
                                   "random_walk",
                                   culture=culture,
                                   params=neuron_params,
                                   num_neurites=1
                                   )
    NetGrowth.Simulate(sim_length)
    # fig, ax = plt.subplots()
    # NetGrowth.plot.PlotNeuron(gid=range(experiment_params["num_neurons"]),
                              # culture = culture, soma_color="k",
                              # axon_color='g', axis=ax, show=True)
    NetGrowth.SaveJson(filepath=save_path)
    NetGrowth.SaveSwc(filepath=save_path, swc_resolution=1)
    # NetGrowth.SaveJson(filepath=tmp_dir)
    # NetGrowth.PlotNeuron(show_nodes=True)
    NetGrowth.ResetKernel()


def Test(neuron_params, sim_length=300, n_samples=10, plot=False):
    swc_folder = os.path.join(os.getcwd(), "tmp_measure")
    CleanFolder(swc_folder)
    RunNetGrowth(n_samples, sim_length, 5,  neuron_params, swc_folder)
    population, _ = NetGrowth.PopulationFromSwc(swc_folder)
    # # return ensembles
    # # rw_corr.plot_results(ensembles, plot=True)
    # info =InfoFromJson(os.path.join(folder,"info.json"))
    # fit = fits[list(fits.keys())[0]]
    # fit = analyse_fit(fit, info=info)
    # fit= OnlyValues(fit)
    # if plot:
    # CleanFolder(folder)
    # RunNetGrowth(1, 1, neuron_params, folder,True)
    # CleanFolder(folder,make=False)
    # print(" ################################### \n")
    # print(" Memory Tau: {} um \n".format(fit["memory"]))
    # print(" Correlation Tau: {} um \n".format(fit["pers_gauss"]))
    # print(" Sigma: {} \n".format(fit["sigma"]))
    # print(" Tortuosity: {} \n".format(fit["tortuosity_local"]))
    # print(" Persistence Lenght: {} um \n".format(fit["pers_length"]))
    # print(" Persistence Lenght from cosine: {} um \n".format(fit["cosine"]))
    # print(" ################################## \n")
    # return fit["pers_length"], fit["tortuosity_local"], fit['cosine']

    # # json.dump(fit,open(os.path.join(folder,"fit.json"),'w'))


def SmoothAngle(sequence, char_length):
    """
    compute the last angle as the weighted average of previous angles,
    the weight is uniform and the characteristic length can be set.

    Parameters:
    ===========

    Sequence to smooth: 1D np.array
    characteristic length: int

    Return:
    ======
    Average angle: float
    """
    average = 0
    for x in range(1, char_length + 1):
        average += sequence[-x]
    return average / (char_length - 1)


def SmoothEnsembleAngle(ensemble, char_length):
    """
    Apply smooth angle for a 2D array, values to average on the second axis
    """
    average_values = []
    for neuron in ensemble:
        last_value = AverageFromHere(neuron.axon.xy[0, :], end_of_the_tube)
        char_length = len(neuron.axon.xy[0, :]) - last_value
        average_values.append(SmoothAngle(neuron.axon.theta, char_length))
    return np.array(average_values)


def AverageFromHere(sequence, x_coord):
    """
    Return the first segment to be considered for averaging the angle.
    Do it checking wheter the x coordinate is greater than the angular cavity
    """
    for n, x in enumerate(sequence):
        if x > x_coord:
            return n
    warnings.warn("last step is reached")
    return n - 1


def ComputeMoment(hist, bin_width):
    tot = sum(hist[0])
    values = hist[0] / tot
    x = hist[1][:-1] + bin_width / 2.
    return sum(x * x * values)


if __name__ == "__main__":

    neuron_params = {
        "filopodia_wall_affinity": 5.,
        "filopodia_finger_length": filopodia_length,
        "filopodia_angular_resolution": 40,
        "axon_angle": 0.,
        "use_tubulin": False,
        "rw_persistence_length": 20.,
        # "rw_delta_corr": 0.1,
        # "rw_memory_tau": 0.7,
        "rw_sensing_angle": 0.15,
        "speed_growth_cone": 1.0,
    }
    names = np.logspace(0, 3, 6)
    pops = []
    for name in names:
        neuron_params["filopodia_wall_affinity"] = name
        pops.append(
            Test(neuron_params, n_samples=n_samples, sim_length=sim_length))
    hists = []
    for pop in pops:
        hists.append(np.histogram(SmoothEnsembleAngle(
            pop, 15), bins=30, range=(-1., 1), density=True))

    bin_width = 2 / 30
    width = bin_width / 6
    moments = {}
    fig, (ax1, ax2) = plt.subplots(2)
    for n, hist in enumerate(hists):
        ax1.bar(hist[1][:-1] * 180 / 3.14 + n * width * 180 / 3.14, hist[0] /
                sum(hist[0]), width=width * 180 / 3.14,
                label="affinity:" + str(names[n])[:5])
        ax2.scatter(names[n],ComputeMoment(hist, bin_width))
        moments[names[n]] = ComputeMoment(hist, bin_width)
    import json
    with open(save_path+".json", 'w') as fp:
        json.dump(moments, fp)


# hist1, alpha=0.2, label ="wall_affinitty: 10.")
# hist2, alpha=0.2, label="wall_affinity : 5.")
# hist3, alpha=0.2, label="wall_affinity : 1.")
# hist4, alpha=0.2, label="wall_affinity : 0.1")
# hist5, alpha=0.2, label="wall_affinity : 0.01")
    ax1.set_xlabel("Angles [degree]")
    ax1.set_ylabel("Probability")
    ax2.set_xlabel("Wall Affinity")
    ax2.set_ylabel("Distribution variance")
    ax2.loglog()
    ax1.legend()
    fig.tight_layout()
    fig.savefig(save_path+".pdf",dpi=300,format="pdf")
    plt.show()
