# VeNoM

## Introduction
*VeNoM* is a study of parameters that may affect performance of an approximate subgraph matching algorithm.
Particularly, it focuses on the breadth and depth of neighborhood considered during similarity computation
of two vertices and creates different instances of the algorithm based on it.

Please cite our paper, if you use our source code.
* "VeNoM: Approximate Subgraph Matching with Enhanced Neighbourhood Structural Information. CODS-COMAD'24"


## Repository Structure
The base folder contains binary for four algorithms, each in a separate folder.
The mapping with respect to the algorithm names in the paper to the corresponding code folder is given below.

* VeNoM-(2,1): two_units
* VeNoM-(3,1): three_units
* VeNoM-(2,2): two_pair
* VeNoM-(1,1): one_unit

## Command
For each algorithm, to execute the binary file run the following command.

```./subgraph file1 file2 file3```

## Parameters
* file1: vertex label file
  - format: `vertex_id vertex_label`
* file2: edge file path
  - format: `vertex_id1 vertex_id2`
* file3: file containing paths of query graphs
  - format: `query_vertex_label_file query_edge_file`

Output is to stdout.
The output of each query graph is further processed and sorted based on their chi-squared value.

## Example
A sample Barabási-Albert graph has been provided in the dataset folder along with some queries and the corresponding query file.

Following commands could be run from the base folder on the sample graph for various algorithms.

* VeNoM-(2,1)
  - ```two_units/subgraph dataset/v1k_l50_labels dataset/v1k_m5_edges dataset/qry_files```
* VeNoM-(3,1): three_units
  - ```three_units/subgraph dataset/v1k_l50_labels dataset/v1k_m5_edges dataset/qry_files```
* VeNoM-(2,2): two_pair
  - ```two_pair/subgraph dataset/v1k_l50_labels dataset/v1k_m5_edges dataset/qry_files```
* VeNoM-(1,1): one_unit
  - ```one_unit/subgraph dataset/v1k_l50_labels dataset/v1k_m5_edges dataset/qry_files```


