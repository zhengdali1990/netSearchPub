# netSearchPub
This program simulate network motifs to find oscillators

makefile will generate executable (default name 'test')

The input has the format of: -t <Topology_File> -f <Interaction_Function> -s <Sampling_Type> -p <parameter_partiation_number> -b <parameter_block_size> -l <lower_boundary_of_parameters> -u<upper_boundary_of_parameters> -n <number_of_parameters> -o <whether_export_timeseries>
Topology_File: text file store all the topologies to be simulated, two examples are provided. Format: 1st line: total number of the nodes. 2nd line: total number of topologies to be simulated. 3rd line to the end: the adjacency matrix of each topology is stored in a line (coloumn wise)
Interaction _Function: protein or proteinMDecay: the former utilize hill coeffecient as interaction term and linear decay, the latter use utilize hill coeffecient as interaction term and Michaelis Menten decay.
Sampling_type: log or linear, represent the method of latin hypercube sampling, log sampling using 10 as base.
parameter partiation number: number of blocks of parameter number.
parameter block size: size of each blocks of parameters (the sum of this number equals the total number of parameters).
lower boundary of parameters: lowest sampling value
higher boundary of parameters: highest sampling value
number of parameters: total number of LHS sampling sive
whether export time series: 0: do not export time series data. 1: export dynaimc time series data.

output: 
a numericalmatrix: 
rows: response to different parameters(only oscillators are shown)
columns: from left to right: parameters, amplitudes, maximum derivative, spikiness(duration of time that the amplitude is larger than average over one period), phase difference, period. 

An example:
-t top3nodes_4edge_2.txt -f protein -s linear -p 5 -b 3 9 9 9 3 -l 0 0 0 1 0 -u 10 1000 1 10 10 -n 10000 -o 0
