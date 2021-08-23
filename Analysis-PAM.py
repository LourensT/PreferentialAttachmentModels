# we analyze the PAM

from PAM import PAM

import numpy as np

import matplotlib.pyplot as plt


# from RGCN volume 2, 8.4.11 en 8.4.12
# we know degree distribution behaves like powerlaw with exponent τ = 3 + δ/m > 2

# tau = 2.5 = 3 + (-1 / 2) = 3 + δ/m
m = 2
d = -1
t = 2000

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
plt.title("Tail Distributions of degrees, tau=2.5 (m = {}, d = {}, t = {})".format(m, d, t))

plt.loglog(base=10)
plt.show()

typicalDist = graph.typicalDistanceDistribution()
# plot distance distribution
plt.scatter(x=typicalDist.keys(), y=typicalDist.values(), color='red')
plt.title("Typical distance distribution, tau=2.5")
plt.show()


# tau = 3.5 = 3 + (1 / 2) = 3 + δ/m
m = 2
d = 1
t = 2000

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
plt.title("Tail Distributions of degrees, tau = 3.5 (m = {}, d = {}, t = {})".format(m, d, t))

plt.loglog(base=10)
plt.show()

typicalDist = graph.typicalDistanceDistribution()
# plot distance distribution
plt.scatter(x=typicalDist.keys(), y=typicalDist.values(), color='red')
plt.title("Typical distance distribution, tau=3.5")
plt.show()