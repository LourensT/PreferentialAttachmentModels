# Generating Preferential Attachment Graphs and analyzing their degree distributions.

References to RGCN I/II refer to Random Graphs and Complex Networks by Remco van der Hofstad, respectively volume 1 and 2. 

Here we generated Preferential Attachment Models (PAMs) as defined in chapter 8 of RGCN volume 1. 

From 8.4.11 en 8.4.12, we know degree distribution behaves like powerlaw with exponent $τ = 3 + δ/m > 2$.

The objective of `Analysis-PAM.py` is to inspect three different definitions for a degree distribution. The conventional definition, the size-biased definition (1.2.2 in RGCN I), and the random friend distribution (the degree distribution of a uniformly selected adjacent vertex to a uniformly selected vertex in the network). See the submodule DegreeDistributions for more details. 

See preview.pdf for the resulting plots. 