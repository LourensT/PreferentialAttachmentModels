import networkx as nx
from scipy import stats
import numpy as np
#import itertools
import matplotlib.pyplot as plt

#import tikzplotlib as plt2tikz

from DegreeDistributions.DegreeDistributions import *

class PAM:

    def __init__(self, n):
        self.G = nx.Graph()

        # generate PAM

    def DegreeDistribrution(self, tail=True):
        return DegreeDistribution(self.G, tail=tail)

    def RandomFriendDegreeDistribution(self, tail=True):
        return RandomFriendDegreeDistribution(self.G, tail=tail)

    def SizeBiasedDegreeDistribution(self, tail=True):
        return SizeBiasedDegreeDistribution(self.G, tail=tail)