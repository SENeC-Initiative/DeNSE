import os
from os.path import join, isfile
import numpy as np

__all__=[
        "SplitSwcFile",
        "ImportSwc"
        ]

def ImportSwc(swc_file, population_to_singles=False):
    """
    Import the swc file as a list of 2 dimensional array.
    First dimension is the neurons space, as many as number of neurons in the file.
    Second dimension is the swc params space, which is headed in the np.array
    Third dimension is the data space.

    Returns
    -------
    a list of bitmorph objects.
    """
    gids = SplitSwcFile(swc_file)
    neurons=[]
    hash_path=swc_file.split(".")[0]

    for gid in range(1,gids+1):
        print("read {}/neuron_{}.swc".format(swc_file,gid))
        neurons.append(np.loadtxt(join(hash_path,"neuron_"+str(gid)+".swc")))
    return neurons, gids

def get_axon_path(neuron, plot=False, neuron_gid=1):
    """
    Only import axon, assume that:
        1) there is only one neurite.
        1) there is no branching in the tree.

    it recognizes the input format for neuron:
    * file path to swc file
    * btmorph NeuronMorphology object
    * np.ndarray
    and converts xy lists to xy and polar coordinates
    """
    if isinstance(neuron,np.ndarray):
        if neuron.shape[1]==2:
            xy= neuron[:,:].transpose()
        elif neuron.shape[1]==6:
            xy= neuron[:,2:4].transpose()
        else:
            raise ValueError("Data type not understood, number neuron parameters: {}".format(neuron.shape[1]))


    elif isfile(neuron) and neuron.endswith(".swc"):
        neuron = np.loadtxt(neuron)
        xy= neuron[:,2:4].transpose()
    else:
        import btmorph2
        if isinstance(neuron, btmorph2.NeuronMorphology):
            if plot:
                neuron.plot_1D()
            xy=btmorph2.get_neuron_path(neuron)[:,5:]
    try:
        angles=_angles_from_xy(xy)
    except:
        raise ValueError("xy neuron data ar ill formed. xy is {}", xy.shape)
    modules=_module_from_xy(xy)
    return xy, np.array([modules,angles])

def get_properties(neuron):
        props = neuron['neurons']['0']['axon_params']
        name= "lp_{} mem_{} var_{}".format(props['rw_delta_corr'],
                                    props['rw_memory_tau'],
                                    props['rw_sensing_angle'])
        return name, props


def _angles_from_xy(path):
    angles=[]
    for n in range(1,len(path[0])):
        deltax=path[0][n]-path[0][n-1]
        deltay=path[1][n]-path[1][n-1]
        rad = np.arctan2(deltay,deltax)
        angles.append(rad)
    return _demodularize(np.array(angles)-angles[0])


def _demodularize(angles):
    shift=0
    demodule=np.zeros((len(angles)))
    for n, theta in enumerate(angles):
        if abs(theta-angles[n-1]) > 3.14:
            shift+=np.sign(angles[n-1])*2*np.pi
        demodule[n]=theta+shift
    return demodule


# def remove_modulus(angle, previous):
    # if angle
def _module_from_xy(path):
    modules=[]
    for n in range(1,len(path[0])):
        deltax=path[0][n]-path[0][n-1]
        deltay=path[1][n]-path[1][n-1]
        module = np.sqrt(deltay**2 + deltax**2)
        modules.append(module)
    return modules


def _lines_to_file(neuron, _file):
    w = open(_file,'w')
    for z in neuron:
        w.write(z)


def SplitSwcFile(input_file):
    """
    write single neurons .swc files from many neurons .swc file
    Btmorph requires this.

    Btmorph read single neuron .swc files
    Netgrowth write many neurons swc file.

    """

    f = open(input_file, 'r')
    if not input_file.endswith("morphology.swc"):
        raise ValueError("SwcFile: morphology.swc expected  instead {} got".format(input_file))
    filename=input_file[:-14]

    if not os.path.exists(filename):
        raise ValueError("NetGrowth simulation folder not found")
        os.makedirs(filename)

    neuron=[]
    gid = 0
    stored_data=False
    for line in f:
        if not line.startswith('#') and line.strip():
            stored_data=True
            # print (line)
            neuron.append(line)
        if stored_data and line.startswith('#'):
            # print( "write", neuron)
            gid+=1
            _lines_to_file(neuron,os.path.join(filename,"neuron_"+str(gid))+".swc")
            neuron =[]
    if neuron:
        gid+=1
        _lines_to_file(neuron,os.path.join(filename,"neuron_"+str(gid))+".swc")
    return gid
