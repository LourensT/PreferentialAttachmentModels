# we analyze the PAM

from PAM import PAM

import numpy as np

import matplotlib.pyplot as plt


m = 2
d = -1
t = 1000

graph = PAM(m, d, t)

distr = graph.DegreeDistribrution(tail=True)
distrBiased = graph.SizeBiasedDegreeDistribution(tail=True)
friendBiased = graph.RandomFriendDegreeDistribution(tail=True)
# plot degree distribution
plt.scatter(x=distr.keys(), y=distr.values(), color='red')
print(max(distr.keys()))
plt.scatter(x=distrBiased.keys(), y=distrBiased.values(), color='green')
print(max(distrBiased.keys()))
plt.scatter(x=friendBiased.keys(), y=friendBiased.values(), color='blue')
print(max(friendBiased.keys()))
plt.legend(["Normal Degree Distribution", "size Biased Degree Distribution", "Friend"])
plt.title("Tail Distributions of degrees, m = {}, d = {}, t = {}".format(m, d, t))
plt.show()
