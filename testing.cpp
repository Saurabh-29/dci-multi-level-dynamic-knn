#include <iostream> 
#include <set>
#include <iterator> 
#include "testing.h"
#include "utils.h"
#include <vector> 
#include <assert.h>
// #include <math.h>
#include <unordered_map>
using namespace std; 


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

void process_data(const double* const data)
{
    cout<<"data is      :";
    int i;
    for(i=0;i<10;i++)
    cout<<"  :"<<data[i];

    const double* loc = &(data[5]);
    cout<<"some selected values :"<<loc[0]<<"\t"<<loc[1]<<endl;
    exit(0);
}



struct Node
{
    int id;
    Node *prev;
    // Node **next;
    // vector < Node > *next;
    set<int> s;
};

struct Cmp{
bool operator() (Val a, Val b)
{   
    if (a.val == b.val)
    {
        return a.global_id<b.global_id;
    }
    return a.val<b.val;
}};
#include <type_traits>


int main() 
{
    Node *node = new Node;
    node->prev= NULL;
    // node->next = NULL;
    node->s.insert(10);
    node->s.insert(20);
    node->s.insert(30);
    node->s.insert(40);
    node->s.insert(50);

    Node *node2 = new Node;
    node2->prev= node;
    
    node2->s.insert(1);
    node2->s.insert(2);
    node2->s.insert(3);
    node2->s.insert(4);

    set<int >::iterator itr1; 
    for (itr1 = node->s.begin(); itr1 != node->s.end(); itr1++)  
    { 
        cout << *itr1<<" "; 
    } 
    cout<<"our results for test \n";

    for (itr1 = node2->s.begin(); itr1 != node2->s.end(); itr1++)  
    { 
        cout << *itr1<<" "; 
    } 
    cout<<"our results for test \n";

    // node->next.pushback(node2);
    // for(int i=0;i<10;i++){ node->next[i] = node2;}
    
    Node *temp = node2; //node->next[0];
    // node->next = new Node* [12];
    for (itr1 = temp->s.begin(); itr1 != temp->s.end(); itr1++)  
    { 
        cout << *itr1<<" "; 
    } 
    cout<<"our results for new range test \n";


    // for(int i=0;i<10;i++){ node->next[i] = node2;}
    // node->next[0] = node2;
    
    // node->next[20] = node2;
   
    cout<<"our results for new range test \n";

    set<int> s1; 
  
    // insert elements in random order 
    s1.insert(40); 
    s1.insert(30); 
    s1.insert(60); 
    s1.insert(20); 
    s1.insert(50); 
    s1.insert(10); 
    s1.insert(1000000);
    set<int >::iterator itr; 
    cout << "\nThe set s1 is : \n"; 
    for (itr = s1.begin(); itr != s1.end(); itr++)  
    { 
        cout << *itr<<" "; 
    } 
    cout << endl; 
    itr = s1.lower_bound(22);
    // itr = s1.find(30);
    cout<<*itr<<"\t"<<*(--itr)<<"\n";
    cout << "s1.lower_bound(40) : \n" 
         << *s1.lower_bound(22) << "\n"
         << *s1.lower_bound(65)
         << endl; 

    // using Cmp = std::integral_constant<decltype(&cmp), &cmp>;
    set<Val, Cmp>  *s_val;

    s_val = new set<Val, Cmp> [12];
    Val var;
    var.val=0.1; var.global_id = 5; var.local_id = 3;
    s_val[0].insert(var);
    var.val=0.2; var.global_id = 6; var.local_id = 0;
    s_val[0].insert(var);
    var.val=0.5; var.global_id = 7; var.local_id = 2;
    s_val[0].insert(var);
    var.val=0.25; var.global_id = 8; var.local_id = 1;
    s_val[0].insert(var);
    var.val=0.25; var.global_id = 10; var.local_id = 1;
    s_val[0].insert(var);

    cout<<"valus in here are :\n";
    set<Val >::iterator itr3; 
    for (itr3 = s_val[0].begin(); itr3 != s_val[0].end(); itr3++)  
    { 
        cout << (*itr3).val<<" "; 
    } 
    Val var2;
    var2.val=0.6;
    var2.global_id = 7;
    itr3 = s_val[0].find(var2);
    if(itr3 == s_val[0].end())
    {cout<<"element not found in find\n";}
    else{
        
        cout<<"find the value "<< (*itr3).val<<"\n ";
    }
    itr3 = s_val[0].lower_bound(var2);
    if(itr3 == s_val[0].end())
    {cout<<"element not found lower_bound\n";}
    else{
        
        cout<<"\nfind the value lb "<< (*itr3).val<<" "<<(*itr3).global_id;
    }
    set<Val >::iterator itr4; 
    itr4 = s_val[0].upper_bound(var2);
    if(itr4 == s_val[0].end())
    {cout<<"element not found upper_bound\n";}
    else{
        
        cout<<"\nfind the value ub "<< (*itr4).val<<" "<<(*itr4).global_id;
    }
    if(itr3==itr4)
    {cout<<"\nthe values match of itr\n";++itr3;}
    cout<<"Values of itr4, itr3 :"<< (*itr4).val<<" "<<(*itr4).global_id<<"\t"<< (*itr3).val<<" "<<(*itr3).global_id;

    var2.val = 0.02;
    itr3 = s_val[0].lower_bound(var2);
    --itr3;
    // --itr3;
    cout<<"\n values at set end lside : "<< (*itr3).val<<" "<<(*itr3).global_id;
    if(itr3 == s_val[0].end())
    {cout<<"\nend of set reached";}
    else{
        
        cout<<"find the value "<< (*itr3).val<<" "<<(*itr3).global_id;
    }



    vector<int> g1; 
    for (int i = 1; i <= 5; i++) 
        g1.push_back(i); 
    
    cout << "\nVector elements are: "; 
    for (auto it = g1.begin(); it != g1.end(); it++) 
        cout << *it << " "; 

    g1[0] = 18;
    cout << "\nVector elements are: "; 
    for (int i = 0; i < 5; i++) 
        cout<<g1[i]<< " ";

    vector<Node *> g2;
    node2->id = 5;
    g2.push_back(node2);
    for (int i = 0; i < 1; i++) 
        cout<<g2[i]->id<< " ";

    node2->id = 7;
    Node *temp2= new Node;
    temp2->id  = node2->id;
    g2.push_back(temp2);
    for (int i = 0; i < 2; i++) 
        cout<<g2[i]->id<< " ";

    node2->id = 6;
    Node *temp3= node2;
    g2.push_back(temp3);

    for (int i = 0; i < 3; i++) 
        cout<<g2[i]->id<< " ";

    cout<<"size of pointers is :"<<sizeof(temp3);
    double* data;
    int dim = 5;
    int intrinsic_dim = 5;
    int num_points = 5;
    int num_queries = 0;
    int num_indices =2;
    assert(posix_memalign((void **)&data, 64, sizeof(double)*dim*(num_points+num_queries)) == 0);

    gen_data(data, dim, intrinsic_dim, num_points+num_queries);
    process_data(data);

    cout << "\n elements are: "; 
    for(int i=0;i<25;i++)
        data[i] = i;
        // cout<<data[i]<<" ";

    double* proj_vec; 
    assert(posix_memalign((void **)&(proj_vec), 64, sizeof(double)*dim*num_indices) == 0);

    cout << "\n pv elements are: "; 
    for(int i=0;i<10;i++)
        proj_vec[i]=1;
        // cout<<proj_vec[i]<<" ";

    double *data_proj;
    assert(posix_memalign((void **)&data_proj, 64, sizeof(double)*num_points*num_indices) == 0);
    matmul(num_points, num_indices, dim, data, proj_vec, data_proj);

    cout << "\n dp elements are: "; 
    for(int i=0;i<10;i++)
        cout<<data_proj[i]<<" ";

    unordered_map<int, int> mymap;
    mymap[2]++;mymap[4]++;
    cout<<"\nvalue of maps "<<mymap[2]<<"\t"<<mymap[4]<<endl;
    return 0; 
}