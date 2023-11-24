#ifndef CONST_H_
 #define CONST_H_

 #define NUM_OF_SYMBOLS 5  // Impacts vertex.* files, stands for s0-s2 symbols
 #define ORDER_CONSTANT 10 // Keeps the factor of maximum size of the internal heap that is considered
 #define MAX_PH_SIZE 20000   // Define maximum size possible for the primary heap
 #define EMPTY_LABEL "__PHI__"  // Only added to the one hop label distribution, does not affect the degree of the vertex
 #define MIN_DEGREE 2   // Impacts vertex.cpp, add dummy neigbour label to nodes
 
 #define TOPK 10  // Defines the highest valued chi squared vertex candidates considered in the greedy approach

 #define LAPLACIAN_BIAS 0.0001 // Corrective laplacian bias for chi square

#endif
