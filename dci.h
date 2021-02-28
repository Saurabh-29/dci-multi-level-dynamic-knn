#ifndef DCI_H
#define DCI_H

#include <iostream> 
#include <set>
#include <iterator> 
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstdio>
#include <vector>
using namespace std; 


typedef struct Val
{
    double val; //to store curr_distance or data_proj, act as key ordering criteria in set as it happens in B tree
    int global_id =0; // store the global id of point 
    int local_id = 0; //stores the local_id of the global_id, currently used in deletion and query
}Val;


struct Cmp{ //Comparator to sort the points on base of val and tie-breaker based on global-id
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


struct Cmp2{ //extra tie-breaker for local_id as global-id can match in case of priority_index when query exactly matches a data point
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

//Node is the data needed for a single point
typedef struct Node
{
    int global_id;      // Global_id of the node
    int level =0;       // Level of the node in tree
    Node *prev;         // Parent of the node   
    vector < Node *> next;      //vector of pointer to the children of node
    int finest_level_points=1;  // Number of points in this node subtree
    const double* loc;          // Stores the pointer to data of the node
    const double* proj_loc;     // Stores the pointer to data_projection
    set<Val, Cmp> top_candidates;       // Stores the top_candidates at level-1 for the node, improves addition/deletion
    set< Val, Cmp>  *s;                 // List(num_indices) of sets to query the data_projection of the children 
}Node;

typedef struct dci {
    int dim;                // (Ambient) dimensionality of data
    int num_comp_indices;   // Number of composite indices
    int num_simp_indices;   // Number of simple indices in each composite index
    int num_points;         // Number of points in the tree
    int num_levels;         // Highest level in the tree
    long long next_point_id;     // ID of the next point to be inserted
    vector < Node *> nodes;      // 
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

// get the next closest projection, given the set, left-right iterator and query_proj
Val dci_next_closest_proj(const set<Val, Cmp> &s, set<Val>::iterator &left, set<Val>::iterator &right, const double query_proj);

// does a knn search on a given set and updates the top-candidates
void dci_query_single_point_single_level(const dci* const dci_inst, Node* point_to_consider_next, int actual_level, int num_neighbours, const double* const query, const double* const query_proj, dci_query_config query_config, set<Val,Cmp >&  top_candidates_new);

// does a knn search on for a single query point using bfs(breadth first search)
int dci_query_single_point(const dci* const dci_inst, int actual_level, int num_neighbours, const double* const query, const double* const query_proj, dci_query_config query_config, set<Val,Cmp >& top_candidates, int if_query, int root_index);

// calculates the points for a node recursively using dfs(depth first search)
void get_finest_level_points(Node* point_to_consider_next);

// updates the top-candidates for a node by call dci_query_single_point with correct arguments
void dci_assign_parent(dci*  dci_inst, const int actual_level, int query_pos, const dci_query_config query_config, int root_index);

// deletes a node from tree, given its id and takes query_config in case re-assignment needs dci_query_single_point
void delete_node(dci* const dci_inst, int id_to_del_next, const dci_query_config query_config);

// reassign the parents of all nodes based on top-candidates of the node
void dci_reassign_all_nodes(dci* const dci_inst, int num_points);

// updates the top-candidates of points_to_reassign recursively using dfs based on root_index, i.e. other tree root index
void dci_recompute_top_candidates(dci* const dci_inst, Node* point_to_reassign, const dci_query_config query_config, int root_index);

// initializes the dci instance
void dci_init(dci* const dci_inst, const int dim, const int num_comp_indices, const int num_simp_indices, const int num_levels);

// generate a dci_tree similar to static version 
void dci_add(dci* const dci_inst, const int dim, const int num_points, const double* const data, const int num_levels, const dci_query_config construction_query_config);

// dynamic insertion function, generates a new tree for new points similar to dci_add and then uses dci_recompute and dci_reassign to merge trees
void dci_subsequent_addition(dci* const dci_inst, const int dim, const int num_points, const double* const data, const int num_levels, const dci_query_config construction_query_config);

// query to get the k nearest neighbour given query, updates the top_candidates
void dci_query(dci* const dci_inst, const int k_nn, const int num_points, const double* const query, const dci_query_config construction_query_config, set<Val,Cmp >*  top_candidates_query);

#endif // DCI_H