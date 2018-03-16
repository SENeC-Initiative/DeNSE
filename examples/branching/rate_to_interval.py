import numpy as np
import matplotlib.pyplot as plt


###########################################
# def rate_discrete_bin(n, t, E, c, A):
    # return A*np.exp(-c*t)*np.power(n, (-E))
###########################################


def process_deltat(E,T,A, perf=False):
    """
    Process simulate interval
    """
    mytime, n, interval = 1,1,0
    ibi, events = [],[]
    end_time = END_TIME
    def Delta_t(n, t, E, T, A):
        """
        Time interval expected for Van Pelt rate
        """
        return np.exp(t/T)*np.power(n, (E))/A
    ###Run the process
    while interval+mytime < end_time:
        interval = int(np.random.exponential(Delta_t(n, mytime, E,T,A)))
        events.append(interval + mytime)
        ibi.append(interval)
        n, mytime= n+1, mytime + interval
    return n, events, ibi

def process_baseline_rate(E,T,A, perf=False):
    mytime, n, last_event = 1,1,0
    ibi, events = [],[]
    end_time = END_TIME
    def rate(n, t, E, T, A):
        return A*np.exp(-t/T) * np.power(n, (-E))
    ###Run the process
    while mytime < end_time:
        if np.random.random() < rate(n, mytime, E,T,A):
            events.append(mytime)
            ibi.append(mytime - last_event)
            [np.random.random() for _ in range(n)] if perf else None
            n, last_event= n+1, mytime
        mytime += resolution
    return n, events, ibi

def process_discretebin(E,_,A, perf=False):
    mytime, n, last_event = 1,1,0
    ibi, events = [],[]
    end_time = END_TIME
    def rate(n, t, E, T, A):
        return A/T * np.power(n, (-E))
    ###Run the process
    while mytime < end_time:
        if np.random.random() < rate(n, mytime, E,END_TIME,A):
            events.append(mytime)
            ibi.append(mytime - last_event)
            [np.random.random() for _ in range(n)] if perf else None
            n, last_event= n+1, mytime
        mytime += 1
    return n, events, ibi

def test_distribution(E,T,A):
    """
    Run distribution, choose which commenting samples
    """
    E = E if isinstance(E,list) else [E]
    for e in E:
        print("computing e = {}".format(e))
        samples[e] = {
                    "baseline" :{"func":process_baseline_rate},
                        # "expected":
                                 # np.power(1+A*T*(1-np.exp(-END_TIME/T)),1./e)},
                    "interval":     {"func":process_deltat},
                    # "discretebin":  {"func":process_discretebin,"expected": END_TIME}
                    }
        for key in samples[e]:
            samples[e][key]["branches"] = []
            samples[e][key]["ibi"] = []
            samples[e][key]["event_time"] = []
            for s in range(ensemble_size):
                n, events, ibi = samples[e][key]["func"](e,T,A)
                samples[e][key]["branches"].append(n)
                samples[e][key]["ibi"].extend(ibi)
                samples[e][key]["event_time"].extend(events)
    return samples

def plot_histograms(samples):
    for num, e in enumerate(E):
        fig, (bx, by, bz) = plt.subplots(3, 1)
        bx.text(0.03, 1.08, str(num)+"."+str('A'),
                horizontalalignment='center',
                verticalalignment='center',
                weight='bold', fontsize=16,
                transform=bx.transAxes)
        by.text(0.03, 1.08, str(num)+"."+str('B'),
                horizontalalignment='center',
                verticalalignment='center',
                weight='bold', fontsize=16,
                transform=by.transAxes)
        bz.text(0.03, 1.08, str(num)+"."+str('C'),
                horizontalalignment='center',
                verticalalignment='center',
                weight='bold', fontsize=16,
                transform=bz.transAxes)
        # fig.suptitle('Branch distribution with parameters:\n A = {}, c= {}, E={} '.format(a,1/c,e), fontsize=17)
        bx.set_title("Number of branching events")
        by.set_title("Inter Branching distribution")
        bz.set_title("Time of Branching distribution")
        bx.set_xlabel("num. branching events (s)")
        bx.set_ylabel("Occurences")
        by.set_ylabel("Occurences")
        bz.set_ylabel("Occurences")
        by.set_xlabel("time (s)")
        bz.set_xlabel("time (s)")

        from matplotlib import cm
        colors = (val for val in cm.jet(np.linspace(0,1,4)))
        for enum, key in enumerate(samples[e]):
            assert_is_list(samples[e][key]["ibi"])
            assert_is_list(samples[e][key]["event_time"])
            # import pdb; pdb.set_trace()  # XXX BREAKPOINT
            c = next(colors)
            # bx.axvline(x=samples[e][key]["expected"], color=c, linestyle='--')
            bx.hist(samples[e][key]["branches"],
                    color=c, bins=55, alpha=0.6, label=key)
            by.hist(samples[e][key]["ibi"],
                    color=c, range=(0, END_TIME/2.), bins=5, alpha=0.6)
            bz.hist(samples[e][key]["event_time"],
                    color=c, range=(0, END_TIME), bins=5, alpha=0.6)
        fig.tight_layout()
        bx.legend()
        fig.savefig("rate_vs_interval_e_{}.pdf".format(e),
                    format='pdf', ppi=300)
    plt.show()

def assert_is_list(lista):
    assert(isinstance(lista,list))
    assert(all(isinstance(item, int) for item in lista))
    # print("list verified")


samples = {}
T = 1000.
A = 0.1
END_TIME = 200
ensemble_size =1000
resolution =1
E = [0.1, 0.3, 0.5, 0.7, 0.99]
# E = [0.5]
samples = test_distribution(E,T,A)
plot_histograms(samples)


# def test_perf():
    # """
    # Test performance of interval prediction
    # vs VanPelt rate
    # """
    # interval = np.ndarray((10, 10))
    # vanpelt = np.ndarray((10, 10))
    # for y in range(10):
        # for num, x in enumerate(np.logspace(1, 4, 10)):
            # # print(x)
            # x = int(x)
            # import time
            # time_0 = time.time()
            # for x in range(x):
                # n, events, ibi = process_vanpelt(c, a*2., 0.7, perf=True)
            # time_vp = time.time() - time_0
            # time_0 = time.time()
            # for x in range(x):
                # n, events, ibi = process_deltat(c, a*2., 0.7)
            # time_interval = time.time() - time_0
            # vanpelt[num, y] = time_vp
            # interval[num, y] = time_interval
    # vanpelt = np.mean(vanpelt, axis=1)
    # interval = np.mean(interval, axis=1)
    # plt.scatter(np.logspace(1, 4, 10), vanpelt/interval, c='b')
    # return vanpelt, interval
