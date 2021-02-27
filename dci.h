#include <iostream> 
#include <set>
#include <iterator> 
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstdio>
#include <vector>
// #include "testing.h"
using namespace std; 

typedef struct Val
{
    double val;
    // double *data_ptr;
    int global_id =0;
    int local_id = 0;
}Val;

struct Cmp{
bool operator() (Val a, Val b)
{   
    if (a.val == b.val)
    {   
        // if(a.global_id==b.global_id){
        //     return a.local_id<b.local_id;
        // }
        
        return a.global_id<b.global_id;
    }
    return a.val<b.val;
}};


struct Cmp2{
bool operator() (Val a, Val b)
{   
    if (a.val == b.val)
    {   
        if(a.global_id==b.global_id){
            return a.local_id<b.local_id;
        }
        
        return a.global_id<b.global_id;
    }
    return a.val<b.val;
}};

typedef struct Node
{
    int global_id;
    int level =0;
    Node *prev;
    // Node **next;
    vector < Node *> next;
    int finest_level_points=1;
    const double* loc;
    set<Val, Cmp> top_candidates;
    set< Val, Cmp>  *s;
}Node;

typedef struct dci {
    int dim;                // (Ambient) dimensionality of data
    int num_comp_indices;   // Number of composite indices
    int num_simp_indices;   // Number of simple indices in each composite index
    int num_points;
    int num_levels;
    long long next_point_id;     // ID of the next point to be inserted
    vector < Node *> nodes;
    // double* data;
    double* proj_vec;       // Assuming column-major layout, matrix of size dim x (num_comp_indices*num_simp_indices)
} dci;

typedef struct dci_query_config {
    bool blind;
    // Querying algorithm terminates whenever we have visited max(num_visited, prop_visited*num_points) points or retrieved max(num_retrieved, prop_retrieved*num_points) points, whichever happens first
    int num_to_visit;
    int num_to_retrieve;
    double prop_to_visit;
    double prop_to_retrieve;
    int field_of_view;
    int min_num_finest_level_points;
} dci_query_config;