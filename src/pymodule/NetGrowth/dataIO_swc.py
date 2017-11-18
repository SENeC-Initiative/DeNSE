import os
import warnings
from os.path import join, isfile
import numpy as np
import matplotlib.pyplot as plt

__all__=[
        "SplitSwcFile",
        "ImportSwc"
        "SwcToSegments",
        "SegmentsToNetgrowth"
        "GetAxonPath"
        "GetProperties"
        "FromPopulation"
        ]

class SWC_ensemble(object):
    """
    Store all the neurons in an unique object.
    Create the ensemble adding a population to it.
    Once a population is ready perform analysis on it.
    """

    def __init__(self, description):
        try:
            self.name, self.description =GetProperties(description)
        except:
            self.name, self.description = description["name"], description["description"]
        #### Store a (N,max_len) matrix, where each neuron maintain it's properties
        self.dendrites={"theta" : [],"r":[],"xy":[]}
        self.axon={"theta":[],"r":[],"xy":[]}

    def add_population(self,neurons):
        for neuron in neurons:
            axon, dendrite = GetPath(neuron=neurons[neuron]['data'])
            self.axon["theta"].append(axon[1])
            self.axon["r"].append(axon[2])
            self.axon["xy"].append(axon[0])
            if dendrite is not None:
                self.dendrites["theta"].append(dendrite[1])
                self.dendrites["r"].append(dendrite[2])
                self.dendrites["xy"].append(dendrite[0])
            ### save shortest shape, alle equivalent is required:

    @classmethod
    def FromPopulation(Ensemble,population):
        ensemble=Ensemble(population['info'])
        ensemble.add_population(population['neurons'])
        return ensemble


def ImportSwc(swc_file):
    """
    Import the "morphology.swc" file and return a list of non equal two-dimensional array.
    First dimension is the neurons space, as many as number of neurons in the file.
    Second dimension is the swc params space, see SWC specifications
    Third dimension is the data space.

    Returns
    -------
    neurons: list np array from .swc file.
    gids: number of this files
    """
    gids = SplitSwcFile(swc_file)
    neurons=[]
    hash_path=swc_file.split(".")[0]
    for gid in range(1,gids+1):
        print("read {}/neuron_{}.swc".format(swc_file,gid))
        neurons.append(np.loadtxt(join(hash_path,"neuron_"+str(gid)+".swc")))
    return neurons, gids

def SplitSwcFile(input_file):
    """
    Write single neurons .swc files from many neurons .swc file
    Netgrowth write multiple neurons swc file, but single neuron files are easier to process.

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

def SwcToSegments(input_file, angle_thresh, length_thresh, element_type=[3]):
    """
    Import a single neuron SWC file, and create a new segments for each branching point.
    Return:
    ======
    Paths: an array of same length (length_thresh) paths (xy)

    """


    segments = segment_from_swc(input_file, element_type)
    # print(segments)
    # plt.plot(segments)
    paths=[]
    for seg in segments:
        if len(seg)>2:
            segment = SegmentToPath(seg)
            paths.extend(SegmentFromAngularThresh(segment,angle_thresh))
    paths =np.array(omologate(paths, length_thresh))
    if paths.shape[0] ==0:
        raise ValueError("Segments list is empty")
    return paths

def SegmentsToNetgrowth(paths, name, info):
    """
    Convert a list of segments to an equivale set of Neuron in NetGrowth format
    NetGrowth_data:
    {"neurons":
        {
        1:{"1":neuronID, "data":swc array}
        2:{"2":neuronID, "data":swc array}

        n:{"3":neuronID, "data":swc array}
        }
    , "info":simulation_parameters}
    """
    NetGrowth_data ={}
    NetGrowth_data["info"]={"name":name,"description":info}
    neurons = {}
    for n, path in enumerate(paths):
        neurons[n]={"gid":n, "data":path.transpose()}
    NetGrowth_data["neurons"]=neurons
    return NetGrowth_data


def GetPath(neuron, plot=False, neuron_gid=1, axon_ID = 2, dendrite_ID=3):
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
        try:
            angles=_angles_from_xy(xy)
        except:
            raise ValueError("xy neuron data ar ill formed. xy is {}", xy.shape)
        modules=_module_from_xy(xy)
        return (xy,modules,angles),None


    elif isfile(neuron) and neuron.endswith(".swc"):
        neuron = np.loadtxt(neuron)
        axon = np.where(neuron[:,1] == axon_ID)[0]
        dendrite = np.where(neuron[:,1] == dendrite_ID)[0]
        axon_xy = neuron[axon,2:4].transpose()
        dendrite_xy = neuron[dendrite,2:4].transpose()
        try:
            angles_axon=_angles_from_xy(axon_xy)
        except:
            raise ValueError("xy neuron data ar ill formed. neuron[axon] is {}, while axon_xy is {}".format(neuron[axon].shape, axon_xy.shape))
        try:
            angles_dendrite=_angles_from_xy(dendrite_xy)
        except:
            warnings.warn("no dendrites in this run")
            return (axon_xy, angles_axon, _module_from_xy(axon_xy)), None

        return (axon_xy,_module_from_xy(axon_xy),angles_axon), (dendrite_xy, _module_from_xy(dendrite_xy),angles_dendrite)

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

def GetProperties(neuron):
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




def  segment_from_swc(input_file, element_type):
    """
    From a single neuron swc file select one element type and cut the tree in a set of
    segments without branching.
    Return them in a list
    """

    f = open(input_file, 'r')
    segments=[[]]
    for line in f:
        if not line.startswith('#') and line.strip():
            if int(line.split()[1]) in element_type:
                sample_number= int(line.split()[0])
                parent_sample= int(line.split()[-1])
                if parent_sample == sample_number-1:
                    segments[-1].append(line)
                else:
                    segments.append([])

    return segments

def SegmentToPath(segment):
    """
    Convert the segment from SWC file into path,
    The enriched information will contain the angle difference between two consecutive pieces, the angle is required to clean the path from sudden curves.
    """
    matrix = np.ndarray((3,len(segment)))
    x_0 = float(segment[0].split()[2])
    y_0 = float(segment[0].split()[3])
    x = float(segment[1].split()[2])
    y = float(segment[1].split()[3])
    theta_0 = np.arctan2(-(y-y_0),-(x-x_0))
    for n, line in enumerate(segment):
        x = float(line.split()[2])
        y = float(line.split()[3])
        theta = np.arctan2(-(y-y_0),-(x-x_0))
        matrix[0][n]= x
        matrix[1][n]= y
        matrix[2][n]= theta - theta_0
        x_0 = x
        y_0 = y
        theta_0 = theta
    plt.plot(matrix[0,:],matrix[1,:])
    # print(matrix[2,matrix[2,:]>1])
    # print(matrix[2,:])
    # plt.scatter(matrix[0,matrix[2,:]>1],matrix[1,matrix[2,:]>1], c='g')
        # if matrix[2][n] > 0.5:
            # return matrix, segment[n:]
    return matrix

def SegmentFromAngularThresh(segment,thresh):
    """
    Cut the segment into subsegments when the angle is bigger
    than threshold
    """
    if np.any(np.abs(segment[2,:])>thresh):
        breakpoints = np.where(np.abs(segment[2,:])>thresh)[0][1:]
        broken=[]
        start=0
        stop=0
        for stop in breakpoints:
            broken.append(segment[:,start:stop-1])
            start=stop
        broken.append(segment[:,stop:])
        return broken
    else:
        return [segment]

def omologate(segments, thresh):
    """
    cut all the segments to 'thresh' length and move each one to start in (0,0)
    remove angle information too
    """
    paths=[]
    for n, segment in enumerate(segments):
        if segment.shape[1]> thresh:
            paths.append((np.subtract(segment.transpose() , segment[:,0]).transpose())[:2,:thresh])
    return paths

