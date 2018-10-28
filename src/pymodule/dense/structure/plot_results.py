import numpy as np
from matplotlib import colors as mcolors
import os

__all__=["PlotRWAnalysis"]

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in colors.items())
sorted_names=[
'b',
'g',
'r',
'c',
'm',
'y',
'k',
'c',
'm',
'y',
'k',
'w'
]

def yerr(my_array):
    return [my_array[:,0]-my_array[:,2],my_array[:,1]-my_array[:,0]]

def PlotRWAnalysis(ensembles, save_path=None, plot=False, axes=None, error_every=10):
    """
    Plot the path, the autocorrelation function, the variance, and the histogram
    for a set of 'experiments'.
    exps expect this format:
    [(path1, 'exp1_name'),...,(path2,'exp2_name')]
    """
    import matplotlib.pyplot as plt
    ev = error_every

    # f, axarr = plt.subplots(2)
    g, ((ax1, ax2), (ax3, ax4))  = plt.subplots(2,2)
    # f,bx=plt.subplots()
    for n, ensemble in enumerate(ensembles):

        # Contraction rate analysis

        ax1.set_title("Contraction rate")
        ax1.set_xlabel(r"length $\mu m$")
        ax1.set_ylabel(r"$\frac{Euclidean}{Length}$|")
        ax1.plot(ensemble.effective_length,
                ensemble.contraction[:,0],
                label=ensemble.name,
                c=sorted_names[n])

        ax1.text(0.03, 1.08,'A',
            horizontalalignment='center',
            verticalalignment='center',
            weight='bold', fontsize = 16,
            transform = ax1.transAxes)

        ax1.errorbar(ensemble.effective_length,
                    ensemble.contraction[:,0],
                    c=sorted_names[n],
                    yerr=yerr(ensemble.contraction),
                    alpha=0.5,
                    errorevery=ev+n, label=ensemble.name)

        ax1.plot(ensemble.effective_length,
                np.exp(- ensemble.effective_length / ensemble.fits['contraction'][0][0])
                * ensemble.fits['contraction'][0][1]+ensemble.fits['contraction'][0][2],
                '--', c='r')

        ax1.legend()

        # MSD analysis

        ax2.set_title(r"\textit{xy} position MSD: $\langle \Delta X^2 \rangle$")
        ax2.set_xlabel(r"length $\mu m$")
        ax2.set_ylabel(r"$\langle X^2(l) \rangle$")
        ax2.plot(ensemble.effective_length,
                ensemble.msd_2D[:,0],
                label=ensemble.name,
                c=sorted_names[n]
                )


        ax2.text(0.03, 1.08,'B',
            horizontalalignment='center',
            verticalalignment='center',
            weight='bold', fontsize = 16,
            transform = ax2.transAxes)
        ax2.errorbar(ensemble.effective_length,
                    ensemble.msd_2D[:,0],
                    yerr=yerr(ensemble.msd_2D),
                    c=sorted_names[n],
                    alpha=0.05,
                    errorevery=ev+n)
        ax2.loglog()

        # if ensemble.msd_2D_fit_wrong == 'msd_2D_quad':
        ballistic =np.where(ensemble.effective_length
                 > ensemble.fits['cosine'][0][1])
        # ##is linear, since the 2 fits are equal
        ax2.plot(ensemble.effective_length[ballistic],
                 ensemble.effective_length[ballistic],
                '--', c='r')
        # else:
            # ## is quadratic and the linear is wrong!
        ax2.plot(ensemble.effective_length  ,
                 ensemble.effective_length**2,
                 '--', c='r')

        ax2.legend()

        # Correlation analysis

        ax3.set_title(r"Correlation function $\langle \cos \theta \rangle$")
        ax3.set_xlabel(r"length $\mu m$")
        ax3.set_ylabel(r"$\langle \bar b_i \bar b_{i+n} \rangle$")
        ax3.plot(ensemble.effective_length,
                ensemble.cosine[:,0],
                label=ensemble.name,
                c=sorted_names[n])
        ax3.plot(ensemble.effective_length  ,
                ensemble.fits['cosine'][0][2]+
                ensemble.fits['cosine'][0][1]*
                np.exp(-ensemble.effective_length/
                ensemble.fits['cosine'][0][0])
                 , '--', c='r')

        ax3.text(0.03, 1.08,'C',
            horizontalalignment='center',
            verticalalignment='center',
            weight='bold', fontsize = 16,
            transform = ax3.transAxes)

        ax3.errorbar(ensemble.effective_length,
                    ensemble.cosine[:,0],
                    c=sorted_names[n],
                    yerr=yerr(ensemble.cosine),
                    alpha=0.05,
                    errorevery=ev+n, label=ensemble.name)

        ax3.semilogy()
        ax3.set_ylim(0.1 , 1)

        ax3.legend()

        # ax3.errorbar(ensemble.effective_length,
                    # ensemble.cosine[:,0],
                    # yerr=yerr(ensemble.cosine),
                    # c=sorted_names[n],
                    # alpha=0.5,
                    # errorevery=10+n)
        # ax3.plot(ensemble.effective_length,np.ones(len(ensemble.effective_length))*
        #     ensemble.fits['tortuosity_local'][0][0], '--', c=sorted_names[n], alpha=0.3)
        # ax3.plot(ensemble.effective_length,\
                  # np.ones(len(ensemble.effective_length))*ensemble.fits['cosine'][0][0], '--', c=sorted_names[n])

        # variance = [path[x]*path[x] for x in np.arange(10, len(path),10)]
        ax4.text(0.03, 1.08,'D',
            horizontalalignment='center',
            verticalalignment='center',
            weight='bold', fontsize = 16,
            transform = ax4.transAxes)
        ax4.set_title(r"Angle MSD: $\langle \Delta\theta^2 \rangle$")
        ax4.set_xlabel(r"length $\mu m$")
        ax4.set_ylabel(r"$\langle \theta^2(l)\rangle$")
        def get_transient(theta):
            for n,x in enumerate(theta):
                if x > 4:
                    break
            return n
        theta_max=get_transient(ensemble.msd_1D[:,0])
        ax4.plot(ensemble.effective_length,
                ensemble.msd_1D[:,0],
                label=ensemble.name,
                c=sorted_names[n])
        ax4.errorbar(ensemble.effective_length,
                ensemble.msd_1D[:,0],
                yerr=yerr(ensemble.msd_1D),
                c=sorted_names[n],
                alpha=0.03,
                errorevery=ev+5*n)
        ax4.plot(ensemble.effective_length[0:theta_max],
                ensemble.effective_length[0:theta_max]*ensemble.fits['msd_1D_ramp'][0][0] +
                ensemble.fits['msd_1D_ramp'][0][1], '--',
                c='r')

        ax4.legend()

    # handles, labels = ax1.get_legend_handles_labels()
    # bx.axis("off")
    # bx.legend(handles, labels, loc='center')


    # ax1.legend(loc='right', shadow=True)
    g.tight_layout()
    # g.tight_layout()
    if save_path is not None:
        # f.savefig(os.getcwd()+"/"+save_path+"_positions_correlation.pdf",dpi=300,format="pdf")
        g.savefig(os.getcwd()+"/"+save_path+"_mean_square_disp__histogram.pdf",dpi=300,format="pdf")
        # f.savefig(os.getcwd()+"/"+save_path+"_legend_.pdf",dpi=300,format="pdf")

    if plot:
        plt.show()
