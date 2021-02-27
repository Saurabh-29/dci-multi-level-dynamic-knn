#include <iostream> 
#include <set>
#include <iterator> 
#include "dci.h"
#include "utils.h"
#include <vector> 
#include <assert.h>
#include <queue>  
#include <unordered_map>
// #include <math.h>
// #include <stdlib.h>

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
    cout<<"value of somethings check :"<< dci_inst->num_comp_indices<< dci_inst->num_simp_indices<<"  "<<node->global_id<<endl;

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
    // cout<<"calling the ----- oooo"<<endl;
    Val var2;
    if((left==s.end()) && (right==s.end())){
        var2.local_id = -1;
        // cout<<"both ends"<<endl;
        // cout<<var2.local_id<<"\t"<<var2.val<<"\t"<<var2.global_id<<endl;
    }else if(left==s.end()){
        var2 = *right;right++;
        // cout<<"left ends"<<endl;
        // cout<<var2.local_id<<"\t"<<var2.val<<"\t"<<var2.global_id<<endl;
    }else if(right==s.end()){
        var2 = *left;left--;
        // cout<<"right ends"<<endl;
        // cout<<var2.local_id<<"\t"<<var2.val<<"\t"<<var2.global_id<<endl;
    }else if((*right).val-query_proj < query_proj -(*left).val){
        var2 = *right;right++;
        // cout<<"no ends r :"<<query_proj<< endl;
        // cout<<var2.local_id<<"\t"<<var2.val<<"\t"<<var2.global_id<<endl;
    }else{
        var2 = *left;left--;
        // cout<<"no ends l :"<<query_proj<< endl;
        // cout<<var2.local_id<<"\t"<<var2.val<<"\t"<<var2.global_id<<endl;
    }
    return var2;
}



void dci_query_single_point_single_level(const dci* const dci_inst, Node* point_to_consider_next, int actual_level, int num_neighbours, const double* const query, const double* const query_proj, dci_query_config query_config, set<Val,Cmp >&  top_candidates_new)
{   cout<<"\n in the final query loop";
    int i, k, m, h;
    double cur_dist;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    int num_candidates = 0;
    int num_points = point_to_consider_next->next.size();
    int num_points_to_retrieve = max_i(query_config.num_to_retrieve, (int)ceil(query_config.prop_to_retrieve*num_points));
    int num_projs_to_visit = max_i(query_config.num_to_visit*dci_inst->num_simp_indices, (int)ceil(query_config.prop_to_visit*num_points*dci_inst->num_simp_indices));
    
    cout<<"\n##search_values   :"<<num_points_to_retrieve<<"\t"<<num_projs_to_visit<<endl;

    double farthest_dists[dci_inst->num_comp_indices];

    for (m = 0; m < dci_inst->num_comp_indices; m++) {
            farthest_dists[m] = 0.0; }
    unordered_map<int, int> counts[num_indices];
    unordered_map<int, int> candidate_dists;

    set<Val, Cmp >::iterator left_pos[num_indices];
    set<Val, Cmp >::iterator right_pos[num_indices];

    Val var;
    var.global_id = 0;
    for(auto it= point_to_consider_next->s[0].begin();it!=point_to_consider_next->s[0].end();it++)
    {
        cout<<(*it).global_id<<"\t"<<(*it).val;
    }
    for(i=0;i<num_indices;i++)
    {   
        var.val = query_proj[i];
        left_pos[i] = point_to_consider_next->s[i].lower_bound(var);
        right_pos[i] = point_to_consider_next->s[i].upper_bound(var);
        if(left_pos[i]==right_pos[i])
            --left_pos[i];
    }
    cout<<"searching for the next-projection"<<endl;
    set<Val, Cmp2> priority_indices[dci_inst->num_comp_indices];
    set<Val, Cmp2 >::iterator prior_itr[dci_inst->num_comp_indices];
    int current_local_pnts[num_indices];
    for (m = 0; m < dci_inst->num_comp_indices; m++){
        for (h = 0; h < dci_inst->num_simp_indices; h++) {
            k = m*dci_inst->num_simp_indices+h;
            // cout<<"inserted_details pre "<<k<<"\t"<<m<<"\t"<<h<<"\t"<<var.val<<"\t"<<var.local_id<<"\t"<<var.global_id<<"\t"<<(*left_pos[k]).val<<"\t"<<(*right_pos[k]).val<<endl;

            var = dci_next_closest_proj( point_to_consider_next->s[k], left_pos[k], right_pos[k], query_proj[k]);
            // cout<<"\n got global_index match in assign :"<<var.global_id;
            current_local_pnts[k]= var.local_id;
            // cout<<"inserted_details "<<k<<"\t"<<m<<"\t"<<h<<"\t"<<var.val<<"\t"<<var.local_id<<"\t"<<var.global_id<<"\t"<<(*left_pos[k]).val<<"\t"<<(*right_pos[k]).val<<endl;
            // cout<<"values abs_d  :"<<var.val<<"\t"<<query_proj[k]<<"\t"<<abs_d(var.val-query_proj[k])<<endl;
            var.val= abs_d(var.val-query_proj[k]);
            var.local_id = h;
            priority_indices[m].insert(var);
            // cout<<" m and h "<<m <<"\t"<<h<<"\t"<<priority_indices[m].size()<<"\t"<<var.global_id<<var.val<<endl;
        }
    }
    
    k = 0;
    // cout<<"\nsearch iterations :" << dci_inst->num_points<<"\t"<<dci_inst->num_simp_indices<<"\t"<<point_to_consider_next->next.size()<<endl;
    while (k < point_to_consider_next->next.size()*dci_inst->num_simp_indices) {
        for (m = 0; m < dci_inst->num_comp_indices; m++) {
            prior_itr[m] = priority_indices[m].begin();
            i = m*dci_inst->num_simp_indices + (*prior_itr[m]).local_id;
            h = (*prior_itr[m]).global_id;
            counts[m][h]++;
            // cout<<"\n got global_index match while search :"<<h<<"\t"<<(*prior_itr[m]).local_id<<"\t"<<counts[m][h]<<"\t"<<m<<"\t"<<i<<endl;
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
            // cout<<"num of candidates are :"<<num_candidates<<"\t"<<priority_indices[m].size()<<endl;
            priority_indices[m].erase(priority_indices[m].begin());
            // cout<<"erased or not??"<<endl;
            var = dci_next_closest_proj( point_to_consider_next->s[i], left_pos[i], right_pos[i], query_proj[i]);

            if(var.local_id==-1){
                 current_local_pnts[i] =-1;
                //  cout<<"----------------------------   "<<priority_indices[m].size()<<endl;
            }
            else{
                current_local_pnts[i]= var.local_id;
                var.val= abs_d(var.val-query_proj[k]);
                var.local_id = i%dci_inst->num_simp_indices;
                priority_indices[m].insert(var);
                // cout<<"size of priority_index :"<<priority_indices[m].size()<<endl;
            }
        }
        if (num_candidates >= num_neighbours)
            {   
                    if (k + 1 >= num_projs_to_visit || num_candidates >= num_points_to_retrieve) {
                        cout<<"\ngoing to break "<<k;
                    break;
                    }
                // return;
            }
        
        k++;
    }
}


void get_finest_level_points(Node* point_to_consider_next){
    int size = point_to_consider_next->next.size();
    for(int i =0; i< size; i++){
        get_finest_level_points(point_to_consider_next->next[i]);
        point_to_consider_next->finest_level_points += point_to_consider_next->next[i]->finest_level_points;
    }   
    return;
}

int dci_query_single_point(const dci* const dci_inst, int actual_level, int num_neighbours, const double* const query, const double* const query_proj, dci_query_config query_config, set<Val,Cmp >& top_candidates, int flag_1=0, int start_node = 0){
    queue< Node * > points_to_expand;
    int query_flag=flag_1;
    Val var;
    cout<<"\n Values before single level query  :";
    for(auto it= top_candidates.begin();it!=top_candidates.end();it++)
    {
        cout<<(*it).global_id<<"\t";
    }
    points_to_expand.push(dci_inst->nodes[start_node]);
     while (!points_to_expand.empty())
     {
        Node* point_to_consider_next = points_to_expand.front();
        cout<<"\n points to expand_next before query  :"<<point_to_consider_next->global_id<<endl;
        //add a logic to handle empty sets in this code
        for(auto it= point_to_consider_next->s[0].begin();it!=point_to_consider_next->s[0].end();it++)
        {
            cout<<(*it).global_id<<"\t";
        }
        cout<<"\n printed... :";
        if(point_to_consider_next->s[0].empty()){
            var.global_id = point_to_consider_next->global_id;
            var.val = compute_dist( dci_inst->nodes[var.global_id+1]->loc, query, dci_inst->dim);
            top_candidates.insert(var);
            points_to_expand.pop();
            cout<<"----------------------------------------------------------------------------------------------------------------------****************************------------------------------------------------------";
        }else{
            set <Val, Cmp>  top_candidate_new;
            dci_query_single_point_single_level(dci_inst, point_to_consider_next, actual_level, num_neighbours, query, query_proj, query_config, top_candidate_new);
            points_to_expand.pop();

            cout<<"query flag is <<"<<query_flag<<endl;
            if(query_flag){
                cout<<"points expanded "<<point_to_consider_next->level<<"\t"<<point_to_consider_next->global_id<<endl;
                if((point_to_consider_next->level-actual_level)>0){
                    for(auto it= top_candidate_new.begin(); it!= top_candidate_new.end(); it++){
                        top_candidates.insert(*it);
                        if((point_to_consider_next->level-actual_level)>1)
                            points_to_expand.push(point_to_consider_next->next[(*it).local_id]);    
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


int dci_assign_parent(dci*  dci_inst, const int actual_level, int query_pos, const dci_query_config query_config) {
    
    int j;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    int cur_num_returned;
    set <Val, Cmp> top_candidate;
    // cout<<"values stored in both-----: "<<&(query[(query_pos)*dci_inst->dim])<<"\t"<<dci_inst->nodes[query_pos+1]->loc<<endl;
    // cur_num_returned = dci_query_single_point(dci_inst, actual_level, 2, &(query[(query_pos)*dci_inst->dim]), &(query_proj[query_pos*num_indices]), query_config, dci_inst->nodes[query_pos+1]->top_candidates);
    // cur_num_returned = dci_query_single_point(dci_inst, actual_level, 2, dci_inst->nodes[query_pos+1]->loc, &(query_proj[query_pos*num_indices]), query_config, dci_inst->nodes[query_pos+1]->top_candidates);
    cur_num_returned = dci_query_single_point(dci_inst, actual_level, 2, dci_inst->nodes[query_pos+1]->loc, dci_inst->nodes[query_pos+1]->proj_loc, query_config, dci_inst->nodes[query_pos+1]->top_candidates);
    cout<<"\n Values after assign tytytt  :";
    for(auto it= dci_inst->nodes[query_pos+1]->top_candidates.begin();it!=dci_inst->nodes[query_pos+1]->top_candidates.end();it++)
    {
        cout<<(*it).global_id<<"\t"<<(*it).val<<"\t";
    }
    assert(cur_num_returned == 1);
    auto it= dci_inst->nodes[query_pos+1]->top_candidates.begin();
    return (*it).global_id;
}

void delete_node(dci* const dci_inst, Node* point_to_delete_next, const dci_query_config query_config){

    int size = point_to_delete_next->next.size();  
    cout<<"\nDeleting node "<< point_to_delete_next->global_id;
    //reassign the parents of the deleted nodes
    //TODO: update finest
    for(int i=0; i<size; i++){
        if(point_to_delete_next->next[i]==NULL)continue;
        cout<<"\n g_id "<<point_to_delete_next->next[i]->global_id<<"\t";
        auto it = point_to_delete_next->next[i]->top_candidates.begin();
        int flag=1;
        ++it;
        for(it; it!=point_to_delete_next->next[i]->top_candidates.end(); it++){
            cout<<(*it).global_id<<"\t"<<(*it).val<<"\t\t";
            if(dci_inst->nodes[(*it).global_id+1]->prev != NULL){
                insert_val_in_par(dci_inst, (*it).global_id+1, point_to_delete_next->next[i]->global_id+1);
                flag=0;
                break;
            }
        }
        
        point_to_delete_next->next[i]->top_candidates.erase(point_to_delete_next->next[i]->top_candidates.begin(), it);
        
        if(flag){
            int par = dci_assign_parent(dci_inst, point_to_delete_next->next[i]->level, point_to_delete_next->next[i]->global_id, query_config);
        }

        it = point_to_delete_next->next[i]->top_candidates.begin();
        for(it; it!=point_to_delete_next->next[i]->top_candidates.end(); it++){
             cout<<(*it).global_id<<"\t"<<(*it).val<<"\t\t";}   
    }
    int global_id = point_to_delete_next->global_id;
    int par = point_to_delete_next->prev->global_id;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    int local_id_in_del = 0;
    auto it = dci_inst->nodes[par+1]->s[0].begin();
    for(int j=0;j< num_indices;j++)
    {
        Val var; 
        var.val = point_to_delete_next->proj_loc[j]; // data_proj[global_id*num_indices + j];
        var.global_id = global_id;
        // var.local_id = size;
        it = dci_inst->nodes[par+1]->s[j].find(var);
        local_id_in_del = (*it).local_id;
        cout<<"\n del data :"<<(*it).global_id<<"\t"<<(*it).val<<"\t"<<(*it).local_id;
        dci_inst->nodes[par+1]->s[j].erase(it);
        // dci_inst->nodes[par]->s[j].insert(var);
    }
    dci_inst->nodes[par+1]->next[local_id_in_del]= NULL;
    dci_inst->nodes[global_id+1]->prev = NULL;
    it = dci_inst->nodes[par+1]->s[0].begin();
    for(it; it!= dci_inst->nodes[par+1]->s[0].end(); it++){
        cout<<"\n remaining data data :"<<(*it).global_id<<"\t"<<(*it).val<<"\t"<<(*it).local_id;
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
    // dci_inst->indices = (btree_p *)malloc(sizeof(btree_p) * num_indices);
    dci_gen_proj_vec(dci_inst->proj_vec, dim, num_indices);
    
    Node* temp = new Node;
    init_node(temp, dci_inst, -1, num_levels);
    dci_inst->nodes.push_back(temp);
    // cout<<"print values <<"<<temp->global_id<<endl;
    // for (i = 0; i < num_indices; i++) {
    //     btree_p_init(&(dci_inst->indices[i]));
    // }
}
void dci_add(dci* const dci_inst, const int dim, const int num_points, const double* const data, const int num_levels, const dci_query_config construction_query_config) {
    int h, i, j;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    double *data_proj;
    assert(posix_memalign((void **)&data_proj, 64, sizeof(double)*num_points*num_indices) == 0);
    //double *data_proj = (double *)memalign(64, sizeof(double)*num_points*num_indices);
    long long next_point_id = dci_inst->next_point_id;
    assert(dim == dci_inst->dim);
    matmul(num_indices, num_points, dci_inst->dim, dci_inst->proj_vec, data, data_proj);
    dci_inst->num_points = num_points;
    Node* temp = new Node[num_points];

    for(j=0;j<num_points;j++)
        {
            init_node(&temp[j], dci_inst, j);
            // cout<<"print values <<"<<(&temp[j])->global_id<<endl;
            dci_inst->nodes.push_back(&temp[j]);
        }
    for(j=0;j<=num_points;j++)
        {
            cout<<"values are :"<<j<<"\t"<<dci_inst->nodes[j]->global_id<<endl;
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
    for(i=0;i<num_points;i++)
    cout<<"level of point "<<i<<"\t"<<data_levels[i]<<endl;

    for(i=0;i< num_levels;i++)
    {
        for(j=0;j< num_points_on_level[i];j++)
        cout<<" point at level "<<i<<" and at "<<j<<"th index is :"<<level_members[i][j]<<endl;
    }

    // confirm the data_projection major, assuming dimnesion of num_indices X num_points

    for(i=0;i<num_points_on_level[actual_num_levels-1];i++)
    {
        j= level_members[actual_num_levels-1][i];
        // dci_inst->nodes[j+1]->level = actual_num_levels-1;
        // dci_inst->nodes[j+1]->prev = dci_inst->nodes[0];
        Val var;
        var.val =0.0; var.global_id =-1; var.local_id = dci_inst->nodes[0]->next.size();
        dci_inst->nodes[j+1]->top_candidates.insert(var);
        dci_inst->nodes[j+1]->loc = &(data[j*dci_inst->dim]);
        dci_inst->nodes[j+1]->proj_loc =  data_proj + j*num_indices;
        insert_val_in_par(dci_inst, 0, j+1);
    }
    cout<<"level n-1 inserted\n"<<endl;

    for(h=actual_num_levels-2;h>=0;h--)
    {
        for(i=0;i<num_points_on_level[h];i++)
        {
            j= level_members[h][i];
            dci_inst->nodes[j+1]->loc =  data + j*dci_inst->dim; // &(data[j*dci_inst->dim]);
            dci_inst->nodes[j+1]->proj_loc =  data_proj + j*num_indices;
            int par =  dci_assign_parent(dci_inst, h, j, construction_query_config);
            // int par = level_members[h+1][rand() % num_points_on_level[h+1]];  //get_parent(j);
            // dci_inst->nodes[j+1]->level = h;
            // dci_inst->nodes[j+1]->prev = dci_inst->nodes[par+1];
            
            insert_val_in_par(dci_inst, par+1, j+1);
            cout<<"\n ***************parent assigned for node "<<j<<"\t"<<par<<"\n\n";
            for(int ii =0;ii<num_points_on_level[h+1];ii++)
            {   
                cout<<level_members[h+1][ii]<<"\t"<<compute_dist(dci_inst->nodes[j+1]->loc , dci_inst->nodes[level_members[h+1][ii]+1]->loc, dci_inst->dim)<<"\t";
            }
            // break;
        }
    }
    get_finest_level_points((dci_inst->nodes[0]));

    for(i=0;i<dci_inst->nodes.size();i++)
    {
        cout<<"stats for node  "<<i<<"\t"<<dci_inst->nodes[i]->global_id<<"\t"<<dci_inst->nodes[i]->level<<"\t"<<dci_inst->nodes[i]->finest_level_points<<endl;
        if(dci_inst->nodes[i]->prev!=NULL)
        {cout<<"previous global id\t"<<dci_inst->nodes[i]->prev->global_id<<"\n";}
        // for(j=0;j<dci_inst->nodes[i]->next.size();j++)
        //     cout<<dci_inst->nodes[i]->next[j]->global_id<<"\t";
        // cout<<endl;
        for (auto itr = dci_inst->nodes[i]->s[0].begin(); itr != dci_inst->nodes[i]->s[0].end(); itr++)  
            { 
                cout << (*itr).global_id<<" "<< (*itr).local_id<<endl; 
            } 
        // for (auto itr = dci_inst->nodes[i]->s[2].begin(); itr != dci_inst->nodes[i]->s[2].end(); itr++)  
        //     { 
        //         cout << (*itr).global_id<<" "<< (*itr).local_id<<endl; 
        //     } 
        
    }

    cout<<"\n\n\nDeleting now.....\n";
    delete_node(dci_inst, dci_inst->nodes[17], construction_query_config);
    for(i=0;i<dci_inst->nodes.size();i++)
    {
        cout<<"stats for node  "<<i<<"\t"<<dci_inst->nodes[i]->global_id<<"\t"<<dci_inst->nodes[i]->level<<"\t"<<dci_inst->nodes[i]->finest_level_points<<endl;
        if(dci_inst->nodes[i]->prev!=NULL)
        {cout<<"previous global id\t"<<dci_inst->nodes[i]->prev->global_id<<"\n";}
        // for(j=0;j<dci_inst->nodes[i]->next.size();j++)
        //     cout<<dci_inst->nodes[i]->next[j]->global_id<<"\t";
        // cout<<endl;
        for (auto itr = dci_inst->nodes[i]->s[0].begin(); itr != dci_inst->nodes[i]->s[0].end(); itr++)  
            { 
                cout << (*itr).global_id<<" "<< (*itr).local_id<<endl; 
            } 
        // for (auto itr = dci_inst->nodes[i]->s[2].begin(); itr != dci_inst->nodes[i]->s[2].end(); itr++)  
        //     { 
        //         cout << (*itr).global_id<<" "<< (*itr).local_id<<endl; 
        //     } 
        
    }
    
    for(i=1;i<dci_inst->nodes.size();i++)
    {   
        cout<<"\n";
        for(j=1;j<dci_inst->nodes.size();j++)
        {
            cout<<compute_dist(dci_inst->nodes[i]->loc, dci_inst->nodes[j]->loc, dci_inst->dim)<<"\t";
        }
    }

    // compute_dist(&(data[h*dci_inst->dim]), query, dci_inst->dim);

}

void dci_query(dci* const dci_inst, const int dim, const int num_points, const double* const query, const double* const data, const int num_levels, const dci_query_config construction_query_config) {
    set <Val, Cmp> top_candidate;
    double *query_proj;
    int num_indices = dci_inst->num_comp_indices*dci_inst->num_simp_indices;
    assert(posix_memalign((void **)&query_proj, 64, sizeof(double)*num_points*num_indices) == 0);
    matmul(num_indices, num_points, dci_inst->dim, dci_inst->proj_vec, data, query_proj);

    int query_pos =0;
        // cur_num_returned = dci_query_single_point(dci_inst, actual_level, 1, &(query[(query_pos)*dci_inst->dim]), &(query_proj[query_pos*num_indices]), query_config, top_candidate, data);

    int num_returned = dci_query_single_point(dci_inst, 0, 2, &(query[(query_pos)*dci_inst->dim]), &(query_proj[query_pos*num_indices]), construction_query_config, top_candidate, 1);
    for(auto it= top_candidate.begin();it!=top_candidate.end();it++)
    {
        cout<<(*it).global_id<<"\t"<<(*it).val<<"\t";
    }
    cout<<"results returned"<<endl;
}


int main() 
{
    int dim = 25;
    int intrinsic_dim = 5;
    int num_points = 25;
    int num_queries = 5;
    // Assuming column-major layout, data is dim x num_points
    
    double* data;
    assert(posix_memalign((void **)&data, 64, sizeof(double)*dim*(num_points+num_queries)) == 0);
    //double* data = (double *)memalign(64, sizeof(double)*dim*(num_points+num_queries));
    gen_data(data, dim, intrinsic_dim, num_points+num_queries);
    // Assuming column-major layout, query is dim x num_queries
    double* query = data + dim*num_points;

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
    cout<<"no error till this point"<<endl;
    dci_add(&dci_inst, dim, num_points, data, num_levels, construction_query_config);

    dci_query(&dci_inst, dim, num_queries,  query, data, num_levels, construction_query_config);
    cout<<"query done!!!!"<<endl;    
    return 0; 
}