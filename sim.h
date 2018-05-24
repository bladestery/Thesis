//
//  sim.h
//  
//
//  Created by Ben Ruktantichoke on 8/9/17.
//
//

#ifndef sim_h
#define sim_h

#include <stdio.h>

#define LENGTH 50
#define WIDTH 50
#define MAX_NODE 100
#define FAT 0.25
#define MIN_HEIGHT 1.2
#define RADIUS 6
#define GROUP 10

struct node {
    int person;
    int blocked;
    int traversed;
    int checked;
    int idx;
    int marked;
    double stability;
    double reachability;
    double x;
    double x_dest;
    double y;
    double y_dest;
    double timer;
    double height;
    double distance[MAX_NODE];
    double distance_group[GROUP];
    unsigned int num_child;
    struct node *child[MAX_NODE];
    unsigned int num_blockers;
    struct node *blockers[MAX_NODE];
    unsigned int num_parent;
    struct node *parent[MAX_NODE];
    struct node *pp[MAX_NODE];
};

struct graph {
    struct node *coordinate[LENGTH][WIDTH];
    struct node AP;
    unsigned int population;
    struct node *people;
    struct node *rr[MAX_NODE];
    unsigned int num_mirrors;
    struct node *mirrors;
};

struct stat {
    double stability;
    int trials;
};

struct graph *generate_graph_rand(unsigned int num, int ap_x, int ap_y, double ap_height, int num_mirrors);
struct graph *generate_graph_unif(int ap_x, int ap_y, double ap_height);
struct graph *generate_graph_group(int ap_x, int ap_y, double ap_height, int group_size);
void find_group(struct graph *graph, int group_size);
struct graph *generate_graph_poisson(unsigned int num, int ap_x, int ap_y, double ap_height, int num_mirrors, int region);
void visualize_graph(struct graph *graph);
void visualize_stability(struct graph *graph);
void visualize_reachability(struct graph *graph);
void destroy_resources(struct graph *graph);

double check_blockage(struct graph *graph);
double update_blockage_d2(struct graph *graph);
double update_blockage(struct graph *graph);

void find_parents(struct graph *graph);
void find_distance(struct graph *graph);
void sort_parent_distance(struct graph *graph);
void sort_parent_height(struct graph *graph);
void update_parents(struct graph *graph, int t);
void update_parents_depth2(struct graph *graph);
double update_parents_group(struct graph *graph, int t);
double update_parents_stable(struct graph *graph, int t);
double update_parents_stable_height(struct graph *graph, int t);

double greedy_matching(struct graph *graph, int t);
double greedy_matching_depth2(struct graph *graph, double rem, int t);
double group_matching(struct graph *graph, int t);
double group_matching_fair(struct graph *graph, int t);
double stable_matching(struct graph *graph, int t);
double stable_matching_fair(struct graph *graph, int t);
double stable_matching_height(struct graph *graph);
double stable_matching_close(struct graph *graph);
double stable_matching_close_height(struct graph *graph);

double update_greedy(struct graph *graph, int t);
double update_greedy_stable(struct graph *graph, int t);
double update_depth2(struct graph *graph);
double update_group(struct graph *graph, double z, int t);
double update_group_stable(struct graph *graph, double z, int t);
double update_group_fair(struct graph *graph, double z, int t);
double update_stable(struct graph *graph, double z, int t);
double update_stable_stable(struct graph *graph, double z, int t);
double update_stable_fair(struct graph *graph, double z, int t);
double update_stable_height(struct graph *graph, double z);
double update_stable_close(struct graph *graph, double z);
double update_stable_close_height(struct graph *graph, double z);

void update_graph(struct graph *graph);
void update_graph_waypoint(struct graph *graph);
void update_graph_waypoint_group(struct graph *graph);
void shift_index(struct graph *graph);
void sort_height_index(struct graph *graph);
void sort_stability(struct graph *graph);
void sort_reachability(struct graph *graph);

void reset_match(struct graph *graph);
void reset(struct graph *graph);

double maximal_matching(struct graph *graph, int t);
double update_perfect(struct graph *graph, int t);
double update_perfect_stable(struct graph *graph, int t);

double get_stability(struct graph *graph);
double calc_stability(struct graph *graph, int t);
double get_jain(struct graph *graph);
double get_stabl(struct graph *graph);
double get_reach(struct graph *graph);

void init_stat(struct stat *stat);
void save_stat(struct graph *graph, struct stat *stat);

#endif /* sim_h */
