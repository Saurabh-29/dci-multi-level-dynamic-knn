#include <iostream> 
#include <set>
#include <iterator> 
#include "dci.h"
#include "utils.h"
#include <assert.h>

int main() 
{
    int dim = 25;
    int intrinsic_dim = 5;
    int num_points = 2500;
    int num_extra_points =1000;
    int num_queries = 5;
    // Assuming column-major layout, data is dim x num_points
    
    double* data;
    assert(posix_memalign((void **)&data, 64, sizeof(double)*dim*(num_points+num_queries+num_extra_points)) == 0);
    //double* data = (double *)memalign(64, sizeof(double)*dim*(num_points+num_queries));
    gen_data(data, dim, intrinsic_dim, num_points+num_queries+num_extra_points);
    // Assuming column-major layout, query is dim x num_queries
    double* query = data + dim*num_points+ dim*num_extra_points;

    // DCI parameters
    int num_comp_indices = 2;
    int num_simp_indices = 5;
    int num_neighbours = 10;
    int num_levels = 3;
    int construction_field_of_view = 10;
    double construction_prop_to_retrieve = 0.002;
    int query_field_of_view = 100;
    double query_prop_to_retrieve = 0.8;

    dci dci_inst;
    dci_init(&dci_inst, dim, num_comp_indices, num_simp_indices, num_levels);
    dci_query_config construction_query_config;
    
    construction_query_config.blind = false;
    construction_query_config.num_to_visit = -1;
    construction_query_config.num_to_retrieve = -1;
    construction_query_config.prop_to_visit = 1.0;
    construction_query_config.prop_to_retrieve = construction_prop_to_retrieve;
    construction_query_config.field_of_view = construction_field_of_view;

    dci_add(&dci_inst, dim, num_points, data, num_levels, construction_query_config);
    dci_subsequent_addition(&dci_inst, dim, num_extra_points, data+num_points*dim, num_levels, construction_query_config);

    delete_node(&dci_inst, 16, construction_query_config);
    delete_node(&dci_inst, 19, construction_query_config);

    dci_query_config query_config;
    query_config.blind = false;
    query_config.num_to_visit = -1;
    query_config.num_to_retrieve = -1;
    query_config.prop_to_visit = 1.0;
    query_config.prop_to_retrieve = query_prop_to_retrieve;
    query_config.field_of_view = query_field_of_view;

    set< Val, Cmp>  *top_query_candidates;
    top_query_candidates = new set<Val, Cmp> [num_queries];
    dci_query(&dci_inst, num_neighbours, num_queries, query, query_config, top_query_candidates); 

    for(int query_pos = 0; query_pos<num_queries; query_pos++){
        cout<<"Results for query "<<query_pos<<"  :: "<<endl;
        for(auto it= top_query_candidates[query_pos].begin();it!=top_query_candidates[query_pos].end();it++){
            cout<<" id :"<<(*it).global_id<<"\t dist :"<<(*it).val<<"\t";
        }
        cout<<endl;
    }
    return 0; 
}