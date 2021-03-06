# we analyze the PAM

from PAM import PAM

from TikzPlotsFromPython.GenerateTikz import GenerateTikz
import os

import numpy as np

import matplotlib.pyplot as plt


# from RGCN volume 2, 8.4.11 en 8.4.12
# we know degree distribution behaves like powerlaw with exponent τ = 3 + δ/m > 2

# tau = 2.5 = 3 + (-1 / 2) = 3 + δ/m
m = 2
d = -1
t = 100000

BASE_FP = os.getcwd() + "\\plots\\"

print(f"degree_distribution_tau{str(3+(d/m)).replace('.', ',')}.tikz")

graph = PAM(m, d, t)
# get distributions
distr = graph.DegreeDistribrution(tail=True)
distrSizeBiased = graph.SizeBiasedDegreeDistribution(tail=True)
distrfriendBiased = graph.RandomFriendDegreeDistribution(tail=True)

# tikz object
plot = GenerateTikz(os.path.join(BASE_FP, f"degree_distribution_tau{str(3+(d/m)).replace('.', ',')}.tikz"), documentation=f"Degree distributions of PAM (t={t}, delta={d}, m={m}), follow powerlaw (tau={str(3+(d/m)).replace('.', ',')}) by RGCN volume 2, 8.4.11 en 8.4.12. ")
plot_log = GenerateTikz(os.path.join(BASE_FP, f"degree_distribution_tau{str(3+(d/m)).replace('.', ',')}_loglog.tikz"), documentation=f"Loglog Degree distributions of PAM (t={t}, delta={d}, m={m}), follow powerlaw (tau={str(3+(d/m)).replace('.', ',')}) by RGCN volume 2, 8.4.11 en 8.4.12. ")
maximum_degrees = max(max(max(distr.keys()), max(distrSizeBiased.keys())), max(distrfriendBiased.keys()))
plot.setConfiguration(0, maximum_degrees, 0, 1, False, False, xlabel="Degrees", ylabel="$P(X > x)$")
plot_log.setConfiguration(0, maximum_degrees, 0, 1, True, True, xlabel="Degrees", ylabel="$P(X > x)$")

# plot degree distribution
plt.scatter(x=distr.keys(), y=distr.values(), color='red')
plot.addSeries(distr, "Degree distribution")
plot_log.addSeries(distr, "Degree distribution")
plt.scatter(x=distrSizeBiased.keys(), y=distrSizeBiased.values(), color='green')
plot.addSeries(distrSizeBiased, "Size biased degree distribution")
plot_log.addSeries(distrSizeBiased, "Size biased degree distribution")
plt.scatter(x=distrfriendBiased.keys(), y=distrfriendBiased.values(), color='blue')
plot.addSeries(distrfriendBiased, "Random friend degree distribution")
plot_log.addSeries(distrfriendBiased, "Random friend degree distribution")
plt.legend(["Normal Degree Distribution", "size Biased Degree Distribution", "Friend"])
plt.title("Tail Distributions of degrees, tau={} (m = {}, d = {}, t = {})".format(str(3+(d/m)).replace('.', ','), m, d, t))
plt.loglog(base=10)
plt.show()

# get typical distance 
typicalDist = graph.typicalDistanceDistribution(sample=2000)
# tikz object
plot_distance = GenerateTikz(os.path.join(BASE_FP, f"typical_distance_tau{str(3+(d/m)).replace('.', ',')}.tikz"), documentation=f"Distribution of typical distance of PAM (t={t}, delta={d}, m={m}). (implies degree distribution has powerlaw distr tau={str(3+(d/m)).replace('.', ',')}")
plot_distance.setConfiguration(0, max(typicalDist.keys()), 0, max(typicalDist.values()), False, False, xlabel="Typical Distance")
# plot distance distribution
plt.scatter(x=typicalDist.keys(), y=typicalDist.values(), color='red')
plot_distance.addSeries(typicalDist, "Typical Distance")
plt.title(f"Typical distance distribution, tau={str(3+(d/m)).replace('.', ',')}")
plt.show()

# tau = 3.5 = 3 + (1 / 2) = 3 + δ/m
m = 2
d = 1
t = 100000

graph = PAM(m, d, t)
# get distributions
distr = graph.DegreeDistribrution(tail=True)
distrSizeBiased = graph.SizeBiasedDegreeDistribution(tail=True)
distrfriendBiased = graph.RandomFriendDegreeDistribution(tail=True)

# tikz object
plot = GenerateTikz(os.path.join(BASE_FP, f"degree_distribution_tau{str(3+(d/m)).replace('.', ',')}.tikz"), documentation=f"Degree distributions of PAM (t={t}, delta={d}, m={m}), follow powerlaw (tau={str(3+(d/m)).replace('.', ',')}) by RGCN volume 2, 8.4.11 en 8.4.12. ")
plot_log = GenerateTikz(os.path.join(BASE_FP, f"degree_distribution_tau{str(3+(d/m)).replace('.', ',')}_loglog.tikz"), documentation=f"Loglog Degree distributions of PAM (t={t}, delta={d}, m={m}), follow powerlaw (tau={str(3+(d/m)).replace('.', ',')}) by RGCN volume 2, 8.4.11 en 8.4.12. ")
maximum_degrees = max(max(max(distr.keys()), max(distrSizeBiased.keys())), max(distrfriendBiased.keys()))
plot.setConfiguration(0, maximum_degrees, 0, 1, False, False, xlabel="Degrees", ylabel="$P(X > x)$")
plot_log.setConfiguration(0, maximum_degrees, 0, 1, True, True, xlabel="Degrees", ylabel="$P(X > x)$")

# plot degree distribution
plt.scatter(x=distr.keys(), y=distr.values(), color='red')
plot.addSeries(distr, "Degree distribution")
plot_log.addSeries(distr, "Degree distribution")
plt.scatter(x=distrSizeBiased.keys(), y=distrSizeBiased.values(), color='green')
plot.addSeries(distrSizeBiased, "Size biased degree distribution")
plot_log.addSeries(distrSizeBiased, "Size biased degree distribution")
plt.scatter(x=distrfriendBiased.keys(), y=distrfriendBiased.values(), color='blue')
plot.addSeries(distrfriendBiased, "Random friend degree distribution")
plot_log.addSeries(distrfriendBiased, "Random friend degree distribution")
plt.legend(["Normal Degree Distribution", "size Biased Degree Distribution", "Friend"])
plt.title("Tail Distributions of degrees, tau={} (m = {}, d = {}, t = {})".format(str(3+(d/m)).replace('.', ','), m, d, t))
plt.loglog(base=10)
plt.show()

# get typical distance 
typicalDist = graph.typicalDistanceDistribution(sample=2000)
# tikz object
plot_distance = GenerateTikz(os.path.join(BASE_FP, f"typical_distance_tau{str(3+(d/m)).replace('.', ',')}.tikz"), documentation=f"Distribution of typical distance of PAM (t={t}, delta={d}, m={m}). (implies degree distribution has powerlaw distr tau={str(3+(d/m)).replace('.', ',')}")
plot_distance.setConfiguration(0, max(typicalDist.keys()), 0, max(typicalDist.values()), False, False, xlabel="Typical Distance")
# plot distance distribution
plt.scatter(x=typicalDist.keys(), y=typicalDist.values(), color='red')
plot_distance.addSeries(typicalDist, "Typical Distance")
plt.title(f"Typical distance distribution, tau={str(3+(d/m))}")
plt.show()