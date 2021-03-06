from datetime import datetime
import networkx as nx
from scipy import stats
import numpy as np
#import itertools
import matplotlib.pyplot as plt
import random

#import tikzplotlib as plt2tikz

from DegreeDistributions.DegreeDistributions import *

class PAM:

    # model parameters
    m = None
    delta_repr = None # this is delta / m

    # graph state
    current_timestep = 0        

    '''
    initializes the PAM

    @param m >=1, number of edges added each timestep
    @param delta >= -m, determines the growth rule
    @param t >= 1, number of timesteps
    @argument recursive=False, there are no benefits to recursive implementation, it is limited by the recursion depth. 
                                but I'm still keeping it in, sue me.

    @constraint: delta >= -m

    '''
    def __init__(self, m, delta, t, recursive=False):
        assert delta >= -m, "parameter delta has to be >= -m"

        self.delta_repr = delta / m
        self.m = m

        # the representation of the PAM in terms of PAM with m = 1.
        self.G_repr = nx.Graph()

        self.degr_distr = {}

        if recursive:
            # generate PAM_m=1 by recursively growing to t timesteps
            self.recursiveGrowth(t)
        else:
            self.iterativeGrowth(t)

    def growFurther(self, t, recursive=False):
        assert t >= self.current_timestep, "timestep given in recursiveGrowth is less than current timestep."

        if recursive:
            self.recursiveGrowth(t)
        else:
            self.iterativeGrowth(t)


    '''
    recursively grow G_repr until it has self.m*t nodes
    '''
    def recursiveGrowth(self, t):
        assert t >= self.current_timestep, "timestep given in recursiveGrowth is less than current timestep."

        if self.G_repr.number_of_nodes() >= self.m*t:
            self.current_timestep = t

            # get the PAM_m,delta from the PAM_1,m/delta
            self._collapseGraph()
            return
        else:
            # add new vertex (first vertex is called "0")
            self.G_repr.add_node(self.G_repr.number_of_nodes())

            # add the appropriate edge with the heuristic  (8.2.1)
            self._addConnection()

            # recursive call
            self.recursiveGrowth(t)

    '''
    iteratively grow G_repr until it has self.m*t nodes
    '''
    def iterativeGrowth(self, t):
        for i in range(self.G_repr.number_of_nodes(), self.m*t):
            # add new vertex (first vertex is called "0")
            self.G_repr.add_node(i)

            # add the appropriate edge with the heuristic  (8.2.1)
            self._addConnection()

            if i%10000 == 0:
                print(f'At node {i} at time {datetime.now().strftime("%H:%M:%S")}')

        self.current_timestep = t
        self._collapseGraph()
        

    def _addConnection(self):
        # the name of the just added vertex
        firstEndpoint = self.G_repr.number_of_nodes()-1

        # self loop probability
        self.degr_distr[firstEndpoint] = 1 + self.delta_repr

        # select randomly the name of another vertex
        secondEndpoint = random.choices(tuple(self.degr_distr.keys()), k = 1, weights = tuple(self.degr_distr.values()))[0]

        # maintain degree distr
        self.degr_distr[secondEndpoint] += 1

        # update self.G
        self.G_repr.add_edge(firstEndpoint, secondEndpoint)


    ''' 
    Explictely calculate the probabilities of attachment for each node. 
    '''
    def _calculateConnectionPMF(self):
        local_t = self.G_repr.number_of_nodes() - 1
        probabilities = np.zeros(self.G_repr.number_of_nodes())

        # self loop probability  (1 + ??)/(local_t(2 + ??) + (1 + ??))
        probabilities[-1] = (1 + self.delta_repr) / (local_t*(2 + self.delta_repr) + (1 + self.delta_repr))

        # other vertices probabilities:  (vertex_degr + ??)/(local_t(2 + ??) + (1 + ??))
        for i in range(len(probabilities)-1):

            # find degree of vertex with a self loop counting double
            adj = self.G_repr.adj[i] 
            vertex_degr = len(adj)+1 if (i in adj.keys()) else len(adj)

            probabilities[i]  = (vertex_degr + self.delta_repr)/(local_t*(2 + self.delta_repr) + (1 + self.delta_repr))

        
        assert abs(sum(probabilities) - 1) < 0.00001, "PMF does not sum to 1, but to {}".format(sum(probabilities))

        return probabilities

    '''
    (PAM_m,delta)_t is defined in terms of(PAM_1,delta/m)_mt.
    Thus to get the (PAM_m,delta)_t, the alternative representation has to be collapsed.
    '''
    def _collapseGraph(self):
        if self.m == 1:
            self.G = self.G_repr
        else:
            # collapses G_repr into G
            self.G = nx.Graph()

            # add t vertices 
            for i in range(self.current_timestep):
                self.G.add_node(i)
            
            for i in range(self.current_timestep): # i denotes the collapsed vertex
                for j in range(self.m): # j iterates through 
                    vertex_index = i*self.m + j

                    # loop through the adjacent vertices
                    for adj_vertex in self.G_repr.adj[vertex_index].keys():

                        # find the name of the collapsed vertex
                        adj_vertex_collapsed = (adj_vertex // self.m)

                        # add edge to the collapsed graph
                        self.G.add_edge(i, adj_vertex_collapsed)

    def draw(self):
        nx.draw_circular(self.G)
        plt.show()

    def DegreeDistribrution(self, tail=True):
        return DegreeDistribution(self.G, tail=tail)

    def RandomFriendDegreeDistribution(self, tail=True):
        return RandomFriendDegreeDistribution(self.G, tail=tail)

    def SizeBiasedDegreeDistribution(self, tail=True):
        return SizeBiasedDegreeDistribution(self.G, tail=tail)

    '''
    Returns distribution of typical distance:
    - the length of the shortest path between two randomly drawn nodes, given that they are connected

    @param sample: The number of randomly drawn 

    '''
    def typicalDistanceDistribution(self, sample=-1):

        all_shortest_paths = []
        if sample == -1:
            #dictionary of dictionaries dict[source][target] = path
            for source, destinations in nx.algorithms.shortest_path(self.G).items():
                for destination, path in destinations.items():
                    all_shortest_paths.append(path)
        else:
            sources = random.choices(list(self.G.nodes), k=sample)
            targets = random.choices(list(self.G.nodes), k=sample)
            for i in range(sample):
                all_shortest_paths.append(nx.algorithms.shortest_path(self.G, sources[i], targets[i]))

        # calculate pmf
        pmf = {}
        numberOfPaths = 0 #if sample > 0, then this will end up being equal to sample
        for path in all_shortest_paths:
            if (len(path)-1) in pmf:
                pmf[len(path)-1] += 1
            else: 
                pmf[len(path)-1] = 1
            numberOfPaths += 1

        print(numberOfPaths)

        #normalize the histogram (paths currently double counted)
        for key in pmf.keys():
            pmf[key] = pmf[key] / numberOfPaths

        assert abs(sum([v for v in pmf.values()]) - 1) < 0.00001, "pmf does not sum to one!!"

        return pmf

    '''
    Returns size of largest connected component (giant components)

    Note: Strictly speaking, we assume the GRG is highly connected, that is, as n -> \inf, 
    liminf of the ( size of the largest connected component / size of network) > 0.
    '''
    def getSizeOfGiantComponent(self):
        # get sorted list of size of all connected components
        component_sizes = [len(c) for c in sorted(nx.connected_components(self.G), key=len, reverse=True)]
        return component_sizes[0]

    '''
    store the networkx graph object
    '''
    def dumpGraph(self, fp='default'):
        import pickle

        if fp == 'default':
            fp = datetime.now().strftime("%H:%M:%S") + ".graph"
        elif not ".graph" in fp:
            fp = fp + ".graph"

        with open(fp, 'wb') as f:
            pickle.dump(self, f)


if __name__=="__main__":
    pam = PAM(2, 0, 1000)
    pam.draw()