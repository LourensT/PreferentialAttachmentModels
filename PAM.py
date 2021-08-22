import networkx as nx
from scipy import stats
import numpy as np
#import itertools
import matplotlib.pyplot as plt

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

    @constraint: delta >= -m

    '''
    def __init__(self, m, delta, t):
        assert delta >= -m, "parameter delta has to be >= -m"

        self.delta_repr = delta / m
        self.m = m

        # the representation of the PAM in terms of PAM with m = 1.
        self.G_repr = nx.Graph()

        # generate PAM_m=1 by recursively growing to t timesteps
        self.recursiveGrowth(t)

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

    def _addConnection(self):
        # the name of the just added vertex
        firstEndpoint = self.G_repr.number_of_nodes()-1

        # select randomly the name of another vertex
        #  TODO: do this sampling by bulk for efficiency
        secondEndpoint = np.random.choice(np.arange(0, self.G_repr.number_of_nodes()), p=self._calculateConnectionPMF())

        self.G_repr.add_edge(firstEndpoint, secondEndpoint)

    def _calculateConnectionPMF(self):
        local_t = self.G_repr.number_of_nodes() - 1
        probabilities = np.zeros(self.G_repr.number_of_nodes())

        # self loop probability  (1 + δ)/(local_t(2 + δ) + (1 + δ))
        probabilities[-1] = (1 + self.delta_repr) / (local_t*(2 + self.delta_repr) + (1 + self.delta_repr))

        # other vertices probabilities:  (vertex_degr + δ)/(local_t(2 + δ) + (1 + δ))
        for i in range(len(probabilities)-1):
            vertex_degr = len(self.G_repr.adj[i])
            vertex_degr = (vertex_degr + self.delta_repr)/(local_t*(2 + self.delta_repr) + (1 + self.delta_repr))
            probabilities[i] = vertex_degr

        
        print(self.G_repr.number_of_nodes())
        print(self.G_repr.number_of_edges())
        print(probabilities)
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
            # collapse graph
            print("CASES FOR m > 1 not yet implented, giving m=1 representation")

            self.G = self.G_repr
            pass


    def draw(self):
        nx.draw(self.G)
        plt.show()

    def DegreeDistribrution(self, tail=True):
        return DegreeDistribution(self.G, tail=tail)

    def RandomFriendDegreeDistribution(self, tail=True):
        return RandomFriendDegreeDistribution(self.G, tail=tail)

    def SizeBiasedDegreeDistribution(self, tail=True):
        return SizeBiasedDegreeDistribution(self.G, tail=tail)

if __name__=="__main__":
    pam = PAM(1, 1, 50)
    pam.draw()