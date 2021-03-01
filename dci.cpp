#include <iostream> 
#include <set>
#include <iterator> 
#include "dci.h"
#include "utils.h"
#include <vector> 
#include <assert.h>
#include <queue>  
#include <unordered_map>

using namespace std; 
static inline double abs_d(double x) {
    return x > 0 ? x : -x;
}

static inline int min_i(int a, int b) {
    return a < b ? a : b;
}

static inline int max_i(int a, int b) {
    return a > b ? a : b;
}

static void dci_gen_proj_vec(double* const proj_vec, const int dim, const int num_indices) {
    int i, j;
    double sq_norm, norm;
    for (i = 0; i < dim*num_indices; i++) {
        proj_vec[i] = rand_normal();
    }
    for (j = 0; j < num_indices; j++) {
        sq_norm = 0.0;
        for (i = 0; i < dim; i++) {
            sq_norm += (proj_vec[i+j*dim] * proj_vec[i+j*dim]);
        }
        norm = sqrt(sq_norm);
        for (i = 0; i < dim; i++) {
            proj_vec[i+j*dim] /= norm;
        }
    }
}

void init_node(Node* node, dci* const dci_inst, int id =-1, int level = 0){
    node->prev=NULL;
    node->s = new set<Val, Cmp> [dci_inst->num_comp_indices*dci_inst->num_simp_indices];
    node->global_id = id;
    node->level = level;
}
void insert_val_in_par(dci* const dci_inst, int par, int child)
{       
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    int size = dci_inst->nodes[par]->next.size();
    // cout<<"size is :"<<size<<endl;
    dci_inst->nodes[par]->next.push_back(dci_inst->nodes[child]);
    dci_inst->nodes[child]->level = dci_inst->nodes[par]->level-1;

    // dci_inst->nodes[j+1]->level = actual_num_levels-1;
    dci_inst->nodes[child]->prev = dci_inst->nodes[par];
    for(int j=0;j< num_indices;j++)
    {
        Val var; 
        var.val = dci_inst->nodes[child]->proj_loc[j]; //   data_proj[(child-1)*num_indices + j];
        var.global_id = child-1;
        var.local_id = size;
        dci_inst->nodes[par]->s[j].insert(var);
    }
}
Val dci_next_closest_proj(const set<Val, Cmp> &s, set<Val>::iterator &left, set<Val>::iterator &right, const double query_proj)
{   
    Val var2;
    if((left==s.end()) && (right==s.end())){
        var2.local_id = -1;
    }else if(left==s.end()){
        var2 = *right;right++;
    }else if(right==s.end()){
        var2 = *left;left--;
    }else if((*right).val-query_proj < query_proj -(*left).val){
        var2 = *right;right++;
    }else{
        var2 = *left;left--;
    }
    return var2;
}



void dci_query_single_point_single_level(const dci* const dci_inst, Node* point_to_consider_next, int actual_level, int num_neighbours, const double* const query, const double* const query_proj, dci_query_config query_config, set<Val,Cmp >&  top_candidates_new)
{   
    int i, k, m, h;
    double cur_dist;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    int num_candidates = 0;
    int num_points = point_to_consider_next->next.size();
    int num_points_to_retrieve = max_i(query_config.num_to_retrieve, (int)ceil(query_config.prop_to_retrieve*num_points));
    int num_projs_to_visit = max_i(query_config.num_to_visit*dci_inst->num_simp_indices, (int)ceil(query_config.prop_to_visit*num_points*dci_inst->num_simp_indices));
    double farthest_dists[dci_inst->num_comp_indices];
    for (m = 0; m < dci_inst->num_comp_indices; m++) {
            farthest_dists[m] = 0.0; }

    unordered_map<int, int> counts[num_indices];
    unordered_map<int, int> candidate_dists;

    set<Val, Cmp >::iterator left_pos[num_indices];
    set<Val, Cmp >::iterator right_pos[num_indices];

    Val var;
    var.global_id = 0;
    for(i=0;i<num_indices;i++)
    {   
        var.val = query_proj[i];
        left_pos[i] = point_to_consider_next->s[i].lower_bound(var);
        right_pos[i] = point_to_consider_next->s[i].upper_bound(var);
        if(left_pos[i]==right_pos[i])
            --left_pos[i];
    }

    set<Val, Cmp2> priority_indices[dci_inst->num_comp_indices];
    set<Val, Cmp2 >::iterator prior_itr[dci_inst->num_comp_indices];
    int current_local_pnts[num_indices];
    for (m = 0; m < dci_inst->num_comp_indices; m++){
        for (h = 0; h < dci_inst->num_simp_indices; h++) {
            k = m*dci_inst->num_simp_indices+h;

            var = dci_next_closest_proj( point_to_consider_next->s[k], left_pos[k], right_pos[k], query_proj[k]);
            current_local_pnts[k]= var.local_id;
            var.val= abs_d(var.val-query_proj[k]);
            var.local_id = h;
            priority_indices[m].insert(var);
        }
    }
    
    k = 0;
    int n_p = point_to_consider_next->s[0].size(); //point_to_consider_next->next.size()
    while (k < n_p*dci_inst->num_simp_indices) {
        for (m = 0; m < dci_inst->num_comp_indices; m++) {
            prior_itr[m] = priority_indices[m].begin();
            i = m*dci_inst->num_simp_indices + (*prior_itr[m]).local_id;
            h = (*prior_itr[m]).global_id;
            counts[m][h]++;
            if(counts[m][h] == dci_inst->num_simp_indices){
                if(candidate_dists.find(h)== candidate_dists.end()){
                    // cur_dist = compute_dist(&(data[h*dci_inst->dim]), query, dci_inst->dim);
                    cur_dist = compute_dist(dci_inst->nodes[h+1]->loc, query, dci_inst->dim);
                    candidate_dists[h] = cur_dist;
                    if (num_candidates < num_neighbours) {
                        var.global_id =h; var.val= cur_dist; var.local_id = current_local_pnts[i];
                        top_candidates_new.insert(var);
                    }else if(cur_dist < (*(--top_candidates_new.end())).val){
                        auto it = --(top_candidates_new.end());
                        top_candidates_new.erase(it);
                        var.global_id =h; var.val= cur_dist; var.local_id = current_local_pnts[i];
                        top_candidates_new.insert(var);
                    }
                    num_candidates++;
                } else {cur_dist = candidate_dists[h];}
                if (cur_dist > farthest_dists[m]) {
                    farthest_dists[m] = cur_dist;
                }
            }
            priority_indices[m].erase(priority_indices[m].begin());
            var = dci_next_closest_proj( point_to_consider_next->s[i], left_pos[i], right_pos[i], query_proj[i]);

            if(var.local_id==-1){
                 current_local_pnts[i] =-1;
            }
            else{
                current_local_pnts[i]= var.local_id;
                var.val= abs_d(var.val-query_proj[k]);
                var.local_id = i%dci_inst->num_simp_indices;
                priority_indices[m].insert(var);
            }
        }
        if (num_candidates >= num_neighbours)
            {   
                    if (k + 1 >= num_projs_to_visit || num_candidates >= num_points_to_retrieve) {
                        // cout<<"\n"<<k<<"\t"<<num_projs_to_visit<<"\t"<<num_candidates<<"\t"<<num_points_to_retrieve;
                        break;
                    }
            }
        
        k++;
    }
    // cout<<"\n Returning no Candidates :"<<num_candidates<<endl;
}


void get_finest_level_points(Node* point_to_consider_next){
    int size = point_to_consider_next->next.size();
    for(int i =0; i< size; i++){
        get_finest_level_points(point_to_consider_next->next[i]);
        point_to_consider_next->finest_level_points += point_to_consider_next->next[i]->finest_level_points;
    }   
    return;
}

int dci_query_single_point(const dci* const dci_inst, int actual_level, int num_neighbours, const double* const query, const double* const query_proj, dci_query_config query_config, set<Val,Cmp >& top_candidates, int if_query=0, int root_index = 0){
    queue< Node * > points_to_expand;
    int query_flag=if_query;
    Val var;

    points_to_expand.push(dci_inst->nodes[root_index]);
    while (!points_to_expand.empty())
     {
        Node* point_to_consider_next = points_to_expand.front();
      
        if(point_to_consider_next->s[0].empty()){
            var.global_id = point_to_consider_next->global_id;
            var.val = compute_dist( dci_inst->nodes[var.global_id+1]->loc, query, dci_inst->dim);
            top_candidates.insert(var);
            points_to_expand.pop();
        }else{
            set <Val, Cmp>  top_candidate_new;
            dci_query_single_point_single_level(dci_inst, point_to_consider_next, actual_level, num_neighbours, query, query_proj, query_config, top_candidate_new);
            points_to_expand.pop();

            if(query_flag){
                if((point_to_consider_next->level-actual_level)>0){
                    for(auto it= top_candidate_new.begin(); it!= top_candidate_new.end(); it++){
                        top_candidates.insert(*it);
                        if((point_to_consider_next->level-actual_level)>1){
                            points_to_expand.push(point_to_consider_next->next[(*it).local_id]); 
                        }   
                    }
                }
            }
            else{
                if((point_to_consider_next->level-actual_level)>2){
                    for(auto it= top_candidate_new.begin(); it!= top_candidate_new.end(); it++){
                        points_to_expand.push(point_to_consider_next->next[(*it).local_id]);
                        }
                }else if((point_to_consider_next->level-actual_level)==2){
                    for(auto it= top_candidate_new.begin(); it!= top_candidate_new.end(); it++){
                            top_candidates.insert(*it);
                        }
                    }   
                }
        
            }  
     } 
    return 1;
}


void dci_assign_parent(dci*  dci_inst, const int actual_level, int query_pos, const dci_query_config query_config, int root_index=0) {
    
    int j;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    int cur_num_returned;
    set <Val, Cmp> top_candidate;
    
    if(actual_level+2> dci_inst->nodes[root_index]->level){
        
       return ;//(*(dci_inst->nodes[query_pos+1]->top_candidates.begin())).global_id;
    }
    cur_num_returned = dci_query_single_point(dci_inst, actual_level, 3, dci_inst->nodes[query_pos+1]->loc, dci_inst->nodes[query_pos+1]->proj_loc, query_config, dci_inst->nodes[query_pos+1]->top_candidates, 0, root_index);

    // auto it= dci_inst->nodes[query_pos+1]->top_candidates.begin();
    return ; //(*it).global_id;
}

void delete_node(dci* const dci_inst, int id_to_del_next, const dci_query_config query_config){

    Node* point_to_delete_next = dci_inst->nodes[id_to_del_next+1];
    int size = point_to_delete_next->next.size();  
    
    //TODO: update finest
    for(int i=0; i<size; i++){
        if(point_to_delete_next->next[i]==NULL)continue;

        auto it = point_to_delete_next->next[i]->top_candidates.begin();
        int flag=1;
        ++it;
        for(it; it!=point_to_delete_next->next[i]->top_candidates.end(); it++){
            // cout<<(*it).global_id<<"\t"<<(*it).val<<"\t\t";
            if(dci_inst->nodes[(*it).global_id+1]->prev != NULL){
                insert_val_in_par(dci_inst, (*it).global_id+1, point_to_delete_next->next[i]->global_id+1);
                flag=0;
                break;
            }
        }
        point_to_delete_next->next[i]->top_candidates.erase(point_to_delete_next->next[i]->top_candidates.begin(), it);
        
        if(flag){
            dci_assign_parent(dci_inst, point_to_delete_next->next[i]->level, point_to_delete_next->next[i]->global_id, query_config);
            int par = (*(dci_inst->nodes[point_to_delete_next->next[i]->global_id+1]->top_candidates.begin())).global_id;
        }

        it = point_to_delete_next->next[i]->top_candidates.begin();
    }
    int global_id = point_to_delete_next->global_id;
    int par = point_to_delete_next->prev->global_id;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    int local_id_in_del = 0;
    auto it = dci_inst->nodes[par+1]->s[0].begin();
    for(int j=0;j< num_indices;j++)
    {
        Val var; 
        var.val = point_to_delete_next->proj_loc[j];
        var.global_id = global_id;
        it = dci_inst->nodes[par+1]->s[j].find(var);
        local_id_in_del = (*it).local_id;
        dci_inst->nodes[par+1]->s[j].erase(it);
    }
    dci_inst->nodes[par+1]->next[local_id_in_del]= NULL;
    dci_inst->nodes[global_id+1]->prev = NULL;
    delete [] dci_inst->nodes[global_id+1]->s;
    // delete dci_inst->nodes[global_id+1];

}

void dci_reassign_all_nodes(dci* const dci_inst, int num_points){
    for(int i=1; i<=num_points; i++){
        if(dci_inst->nodes[i]->prev==NULL)continue;//this node is deleted but still in vector

        // cout<<"\nCurrent Node : top_candidate : level : current_top   ::  "<< dci_inst->nodes[i]->global_id<< "  : "<< (*(dci_inst->nodes[i]->top_candidates.begin())).global_id<< "  : "<<dci_inst->nodes[i]->level<< "  : "<<dci_inst->nodes[i]->prev->global_id;
        
        if(dci_inst->nodes[i]->prev->global_id == (*(dci_inst->nodes[i]->top_candidates.begin())).global_id){
            // cout<<"  NO REASSIGNMENT REQUIRED";
        }else{
            int new_par = (*(dci_inst->nodes[i]->top_candidates.begin())).global_id;
            int curr_par = dci_inst->nodes[i]->prev->global_id;
            //delete from curr_par and add in new parent
            insert_val_in_par(dci_inst, new_par+1, i);

            //for deletion (maybe write a function to do this?)
            int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
            int local_id_in_del = 0;
            auto it = dci_inst->nodes[curr_par+1]->s[0].begin();

            for(int j=0;j< num_indices;j++){
                Val var; var.val = dci_inst->nodes[i]->proj_loc[j]; var.global_id = dci_inst->nodes[i]->global_id;
                it = dci_inst->nodes[curr_par+1]->s[j].find(var);
                local_id_in_del = (*it).local_id;
                dci_inst->nodes[curr_par+1]->s[j].erase(it);
            }
            dci_inst->nodes[curr_par+1]->next[local_id_in_del]= NULL;//is needed? maybe remove ths in future
        }
    
    }
}

void dci_recompute_top_candidates(dci* const dci_inst, Node* point_to_reassign, const dci_query_config query_config, int root_index)
{
    int size = point_to_reassign->next.size();
    auto it = point_to_reassign->s[0].begin();
    for(it; it!= point_to_reassign->s[0].end(); it++){
        int query_pos = (*it).global_id;
        dci_assign_parent(dci_inst, dci_inst->nodes[query_pos+1]->level, query_pos, query_config, root_index);
        dci_recompute_top_candidates(dci_inst, dci_inst->nodes[query_pos+1], query_config, root_index);
    }
}

void dci_init(dci* const dci_inst, const int dim, const int num_comp_indices, const int num_simp_indices, const int num_levels) {
    
    int i;
    int num_indices = num_comp_indices*num_simp_indices;
    srand48(time(NULL));
    
    dci_inst->dim = dim;
    dci_inst->num_comp_indices = num_comp_indices;
    dci_inst->num_simp_indices = num_simp_indices;
    dci_inst->num_points = 0;
    dci_inst->next_point_id = 0LL;
    dci_inst->num_levels = num_levels;
    
    assert(posix_memalign((void **)&(dci_inst->proj_vec), 64, sizeof(double)*dim*num_indices) == 0);
    dci_gen_proj_vec(dci_inst->proj_vec, dim, num_indices);
    
    Node* points_to_add = new Node;
    init_node(points_to_add, dci_inst, -1, num_levels);
    dci_inst->nodes.push_back(points_to_add);
}

void dci_add(dci* const dci_inst, const int dim, const int num_points, const double* const data, const int num_levels, const dci_query_config construction_query_config) {
    int h, i, j;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    double *data_proj;
    assert(posix_memalign((void **)&data_proj, 64, sizeof(double)*num_points*num_indices) == 0);
    long long next_point_id = dci_inst->next_point_id;

    assert(dim == dci_inst->dim);
    matmul(num_indices, num_points, dci_inst->dim, dci_inst->proj_vec, data, data_proj);
    dci_inst->num_points = num_points;

    // add num_points to the vector of nodes
    Node* points_to_add = new Node[num_points];
    for(j=0;j<num_points;j++){
            init_node(&points_to_add[j], dci_inst, j);
            dci_inst->nodes.push_back(&points_to_add[j]);
        }

    int *data_levels;
    double promotion_prob;
    int num_points_on_level[num_levels];
    int level_relabelling[num_levels];
    data_levels = (int *)malloc(sizeof(int)*num_points);
    promotion_prob = pow((double)num_points, -1.0 / num_levels);
    
    for (i = 0; i < num_levels; i++) {
        num_points_on_level[i] = 0;
    }
    for (j = 0; j < num_points; j++) {
        for (i = 0; i < num_levels - 1; i++) {
            if (drand48() > promotion_prob) {
                break;
            }
        }
        num_points_on_level[i]++;
        data_levels[j] = i;
        }
    h = 0;
    for (i = 0; i < num_levels; i++) {
        if (num_points_on_level[i] > 0) {
            level_relabelling[i] = h;
            h++;
        } else {
            level_relabelling[i] = -1;
        }
    }
    int actual_num_levels = h;
    for (i = 0; i < num_levels; i++) {
            if (level_relabelling[i] >= 0) {
                num_points_on_level[level_relabelling[i]] = num_points_on_level[i];
            }
        }
    int **level_members;
    level_members = (int **)malloc(sizeof(int*)*actual_num_levels);
    for (i = 0; i < actual_num_levels; i++) {
        level_members[i] = (int *)malloc(sizeof(int)*num_points_on_level[i]);
        h = 0;
        for (j = 0; j < num_points; j++) {
            if (level_relabelling[data_levels[j]] == i) {   
                level_members[i][h] = j;
                h++;
            }
        }
        assert(h == num_points_on_level[i]);
    }
  
    // assign the parent for the top-level points to a root-index for ease of access
    for(i=0;i<num_points_on_level[actual_num_levels-1];i++){
        j= level_members[actual_num_levels-1][i];
        Val var;
        var.val =0.0; var.global_id =-1; var.local_id = dci_inst->nodes[0]->next.size();
        dci_inst->nodes[j+1]->top_candidates.insert(var);
        dci_inst->nodes[j+1]->loc = &(data[j*dci_inst->dim]);
        dci_inst->nodes[j+1]->proj_loc =  data_proj + j*num_indices;
        insert_val_in_par(dci_inst, 0, j+1);
    }

    // assign parents for the lower levels
    for(h=actual_num_levels-2;h>=0;h--){
        for(i=0;i<num_points_on_level[h];i++){
            j= level_members[h][i];
            dci_inst->nodes[j+1]->loc =  data + j*dci_inst->dim; // &(data[j*dci_inst->dim]);
            dci_inst->nodes[j+1]->proj_loc =  data_proj + j*num_indices;
            dci_assign_parent(dci_inst, h, j, construction_query_config);
            int par = (*(dci_inst->nodes[j+1]->top_candidates.begin())).global_id;
            
            insert_val_in_par(dci_inst, par+1, j+1);
        }
    }
    // TODO: make better use of finest level points, not used currently
    get_finest_level_points((dci_inst->nodes[0]));

    // free(data_proj);
    // delete [] points_to_add;
}



void dci_subsequent_addition(dci* const dci_inst, const int dim, const int num_points, const double* const data, const int num_levels, const dci_query_config construction_query_config) {
    int h, i, j;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    double *data_proj;
    assert(posix_memalign((void **)&data_proj, 64, sizeof(double)*num_points*num_indices) == 0);
    long long next_point_id = dci_inst->next_point_id;

    assert(dim == dci_inst->dim);
    matmul(num_indices, num_points, dci_inst->dim, dci_inst->proj_vec, data, data_proj);
   
    Node* points_to_add = new Node[num_points+1]; 
    for(j=0;j<num_points;j++){
            init_node(&points_to_add[j], dci_inst, j+dci_inst->num_points);
            dci_inst->nodes.push_back(&points_to_add[j]);
        }

    init_node(&points_to_add[num_points], dci_inst, dci_inst->num_points+num_points, num_levels);
    dci_inst->nodes.push_back(&points_to_add[num_points]);

    int *data_levels;
    double promotion_prob;
    int num_points_on_level[num_levels];
    int level_relabelling[num_levels];
    data_levels = (int *)malloc(sizeof(int)*num_points);
    promotion_prob = pow((double)num_points + (double)dci_inst->num_points/3, -1.0 / num_levels);
    
    for (i = 0; i < num_levels; i++) {
        num_points_on_level[i] = 0;
    }
    for (j = 0; j < num_points; j++) {
        for (i = 0; i < num_levels - 1; i++) {
            if (drand48() > promotion_prob) {
                break;
            }
        }
        num_points_on_level[i]++;
        data_levels[j] = i;
        }
    h = 0;
    for (i = 0; i < num_levels; i++) {
        if (num_points_on_level[i] > 0) {
            level_relabelling[i] = h;
            h++;
        } else {
            level_relabelling[i] = -1;
        }
    }
    int actual_num_levels = h;
    for (i = 0; i < num_levels; i++) {
            if (level_relabelling[i] >= 0) {
                num_points_on_level[level_relabelling[i]] = num_points_on_level[i];
            }
        }
    int **level_members;
    level_members = (int **)malloc(sizeof(int*)*actual_num_levels);
    for (i = 0; i < actual_num_levels; i++) {
        level_members[i] = (int *)malloc(sizeof(int)*num_points_on_level[i]);
        h = 0;
        for (j = 0; j < num_points; j++) {
            if (level_relabelling[data_levels[j]] == i) {   
                level_members[i][h] = j;
                h++;
            }
        }
        assert(h == num_points_on_level[i]);
    }

    int root_index = num_points + dci_inst->num_points+1;
    dci_inst->nodes[root_index]->level = actual_num_levels;
   
    for(i=0;i<num_points_on_level[actual_num_levels-1];i++)
    {
        j= level_members[actual_num_levels-1][i];
        Val var;
        var.val = 10000000.0; var.global_id =root_index-1; var.local_id = dci_inst->nodes[root_index]->next.size();
        dci_inst->nodes[dci_inst->num_points+j+1]->top_candidates.insert(var);
        dci_inst->nodes[dci_inst->num_points+j+1]->loc = &(data[j*dci_inst->dim]);
        dci_inst->nodes[dci_inst->num_points+j+1]->proj_loc =  data_proj + j*num_indices;
        insert_val_in_par(dci_inst, root_index, dci_inst->num_points+j+1);
    }
    
    for(h=actual_num_levels-2;h>=0;h--)
    {
        for(i=0;i<num_points_on_level[h];i++)
        {
            j= level_members[h][i];
            dci_inst->nodes[dci_inst->num_points+j+1]->loc =  data+j*dci_inst->dim; // &(data[j*dci_inst->dim]);
            dci_inst->nodes[dci_inst->num_points+j+1]->proj_loc =  data_proj + j*num_indices;
            dci_assign_parent(dci_inst, h, dci_inst->num_points+j, construction_query_config, root_index);
            int par = (*(dci_inst->nodes[dci_inst->num_points+j+1]->top_candidates.begin())).global_id;
            insert_val_in_par(dci_inst, par+1, dci_inst->num_points+j+1);
        }
    }

    // updates top-candidate of both trees with unseen points
    dci_recompute_top_candidates(dci_inst, dci_inst->nodes[root_index], construction_query_config, 0);
    dci_recompute_top_candidates(dci_inst, dci_inst->nodes[0], construction_query_config, root_index);

    //if actual level is same: add original tree root as top-candidate in level n-1 in new tree
    if(actual_num_levels == dci_inst->num_levels){
        for(i=0;i<num_points_on_level[actual_num_levels-1];i++){
            j= level_members[actual_num_levels-1][i];
            Val var;
            var.val = 0.0; var.global_id =-1; var.local_id = dci_inst->nodes[0]->next.size();
            dci_inst->nodes[dci_inst->num_points+j+1]->top_candidates.insert(var);
        }
    }

    //reassign the parents of all node based on top-candidate value
    dci_reassign_all_nodes(dci_inst, dci_inst->num_points+num_points);
    dci_inst->num_points+= num_points;
    dci_inst->nodes.pop_back();

    // delete [] points_to_add;
    return;
}

void dci_query(dci* const dci_inst, const int k_nn, const int num_queries, const double* const query, const dci_query_config construction_query_config, set<Val,Cmp >*  top_candidates_query) {
    
    double *query_proj;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    assert(posix_memalign((void **)&query_proj, 64, sizeof(double)*num_queries*num_indices) == 0);
    matmul(num_indices, num_queries, dci_inst->dim, dci_inst->proj_vec, query, query_proj);

    int query_pos =0; int actual_level =0; int num_neighbours = k_nn;
    for(query_pos = 0; query_pos<num_queries; query_pos++){
        int num_returned = dci_query_single_point(dci_inst, actual_level, num_neighbours, &(query[(query_pos)*dci_inst->dim]), &(query_proj[query_pos*num_indices]), construction_query_config, top_candidates_query[query_pos], 1);
        auto it= top_candidates_query[query_pos].begin();
        // for(it; it!=top_candidates_query[query_pos].end();it++)
        // {
        //     cout<<(*it).global_id<<"\t"<<(*it).val<<"\t";
        // }
        // it= top_candidates_query[query_pos].begin();
        advance(it, k_nn);
        top_candidates_query[query_pos].erase(it, top_candidates_query[query_pos].end());
    }

    free(query_proj);
}

