//
//  main.c
//  
//
//  Created by Ben Ruktantichoke on 8/9/17.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "sim.h"

#define TRIALS 10000

int main(int argc, char *argv[]) {
    if (argc != 9) {
        fprintf(stderr, "./main %%ui[Number of Nodes] \
                                %%ui[X Dimension] \
                                %%ui[Y Dimension] \
                                %%ui[AP x-coordinate] \
                                %%ui[AP y-coordinate] \
                                %%f[AP Height] \
                                %%ui[Time Step] \
                                %%ui[Num AP - Not Implemented]\n");
        return 1;
    }
    
    unsigned int num = atoi(argv[1]);
    if (num == 0) {
        fprintf(stderr, "Number of nodes must not be 0\n");
        return 1;
    }
    int length = atoi(argv[3]);
    int width = atoi(argv[2]);
    if (length < 5 || width < 5) {
        fprintf(stderr, "width and length must be positive integers > 5\n");
        return 1;
    }
    int ap_x = atoi(argv[4]);
    int ap_y = atoi(argv[5]);
    double ap_height = atof(argv[6]);
    int timestep = atoi(argv[7]);
    if (timestep > 200) {
        fprintf(stderr, "Number of timesteps must be 200 or less\n");
        return 1;
    }
    int num_mirrors = atoi(argv[8]);
    int region = 4;
    time_t t;
    srand((unsigned) time(&t));
    
    int group_size[6] = {1, 2, 4, 5, 8, 10};
    for (int i = 0; i < 6; i++) {
        if (group_size[i] > num) {
            break;
        }
        fprintf(stdout, "Group Size:\n%d\n", group_size[i]);
        double capacity[MAXT] = {0};
        double delay[MAXT] = {0};
        double stability[MAXT] = {0};
        double fair[MAXT] = {0};
        double reach[MAXT] = {0};
        double count = 0;
        double count2 = 0;
        double through_all[MAX_NODE] = {0};
        double latency_all[MAX_NODE] = {0};
        double through[MAX_NODE] = {0};
        double latency[MAX_NODE] = {0};
        double through_max[MAX_NODE] = {0};
        double latency_max[MAX_NODE] = {0};
        for (int x = 0; x < TRIALS; x++) {
            int fail = 0;
            struct node *sorth[MAX_NODE] = {0};
            double temp_ta[MAX_NODE] = {0};
            double temp_la[MAX_NODE] = {0};
            //fprintf(stderr, "x:%d\n", x);
            
            //------------GENERATE GRAPHS--------------
            /*
            struct graph *graph = NULL;
            if (num > 4) {
                graph = generate_graph_poisson(width, length, ap_x, ap_y, num, ap_height, 4);
            } else {
                graph = generate_graph_poisson(width, length, ap_x, ap_y, num, ap_height, 2);
            }
             */
            
            //struct graph *graph = generate_graph_unif(width, length, ap_x, ap_y, num, ap_height);
            struct graph *graph = generate_graph_group(width, length, ap_x, ap_y, ap_height, num, group_size[i]);
            fill_group(graph, group_size[i]);
            sort_group_capacity(graph);
            //sort_group_distance(graph);
            if (graph == NULL) {
                return 1;
            }
            sort_pointer(graph, sorth);
            
            //-------------Check Blockage----------------
            //sort_height_index(graph);
            double ret = check_blockage(graph);
            
            //visualize_graph(width, length, graph);
            find_parents(graph);
            find_distance(graph);
            sort_parent_capacity(graph);
            //sort_parent_distance(graph);
            //sort_parent_height(graph);
            
            //--------------Matching--------------------
            
            //ret = greedy_matching(graph, 0);
            //ret = greedy_matching_depth2(graph, ret, 0);

            //ret = maximal_matching(graph, 0);

            ret = group_matching(graph, 0);
            
            //ret = stable_matching(graph, 0);
            
            //------------Update-Statistics--------------
            
            update_capacity_delay(graph, 0);
            capacity[0] += get_capacity(graph);
            delay[0] += get_delay(graph);
            for (int y = 0; y < num; y++) {
                temp_ta[y] += sorth[y]->capacity;
                temp_la[y] += sorth[y]->delay;
            }

            if (((int)(ret + 0.5)) > 0) {
                fail = 1;
            }
            fair[0] += get_stabl(graph);
            if (fail == 1) {
                reach[0] += get_reach(graph);
            }
            //-----------Simulate Node Mobility-----------
            for (int y = 1; y < timestep; y++) {
                //-----------Update-Graph------------------
                update_graph_waypoint_group(width, length, graph); //random waypoint group
                //shift_index(graph);
                //sort_stability(graph);
                sort_reachability(graph);

                ret = update_blockage(graph);
                sort_group_parent_capacity(graph);

                //update_parents_depth2(graph);
                //update_parents(graph, y);
                //ret = update_parents_stable(graph, y);
                ret = update_parents_group(graph, y);
             
                //reset(graph);
                //find_parents(graph);
                find_distance(graph);
                sort_parent_capacity(graph);
                //sort_parent_distance(graph);
                
                //---------Fair&Stable-Matching-Algorithms------------
                //ret = update_greedy(graph, y);
                //ret = update_depth2(graph);
                //ret = update_greedy_stable(graph, y);

                //ret = update_perfect(graph, y);
                //ret = update_perfect_stable(graph, y);

                //ret = update_stable(graph, ret, y);
                //ret = update_stable_stable(graph, ret, y);
                //ret = update_stable_fair(graph, ret, y);
                
                //ret = update_group(graph, ret, y);
                //ret = update_group_stable(graph, ret, y);
                ret = update_group_fair(graph, ret, y);

                
                //------------Baseline-Matching-Algorithms-------------
                //ret = greedy_matching(graph, y);
                
                //ret = greedy_matching_depth2(graph, ret, y);
             
                //ret = maximal_matching(graph, y);

                //ret = group_matching_fair(graph, y);
                //ret = group_matching(graph , y);
                
                //ret = stable_matching(graph, y);
                //ret = stable_matching_fair(graph, y);
                
                //----------Update-Statistics-----------
                if (((int)(ret + 0.5)) > 0) {
                    fail = 1;
                }
                update_capacity_delay(graph, y);
                capacity[y] += get_capacity(graph);
                delay[y] += get_delay(graph);
                for (int z = 0; z < num; z++) {
                    temp_ta[z] += sorth[z]->capacity;
                    temp_la[z] += sorth[z]->delay;
                }
                
                stability[y] += calc_stability(graph, y);
                fair[y] += get_stabl(graph);
                if (fail == 1) {
                    double tt = get_reach(graph);/*
                    if (tt > 10) {
                        visualize_reachability(width, length, graph);
                        for (int b = 0; b < graph->population; b++) {
                            if (graph->people[b].person == 1) {
                                fprintf(stderr, "x: %d, y: %d, x_dest: %d, y_dest: %d, timer: %d\n", (int) (graph->people[b].x + 0.5),
                                        (int) (graph->people[b].y + 0.5), (int) (graph->people[b].x_dest + 0.5), (int) (graph->people[b].y_dest + 0.5), (int) (graph->people[b].timer + 0.5));
                            }
                        }
                    }*/
                    reach[y] += tt;
                }
            }
            //------Update-throughput-latency-all-trials-----------
            for (int y = 0; y < num; y++) {
                through_all[y] += temp_ta[y]/timestep;
                latency_all[y] += temp_la[y]/timestep;
            }
            
            //--------Count-for-trials-with-unreachable-nodes-------
            if (fail == 1) {
                count++;
                for (int y = 0; y < num; y++) {
                    through[y] += temp_ta[y]/timestep;
                    latency[y] += temp_la[y]/timestep;
                }
            }
            
            //-------Ret-for-trials-with-most-unreachability-------
            ret = get_reach(graph);
            if (ret > count2) {
                count2 = ret;
                for (int y = 0; y < num; y++) {
                    through_max[y] = temp_ta[y]/timestep;
                    latency_max[y] = temp_la[y]/timestep;
                }
            }
            
            destroy_resources(graph);
        }
        
        //---------Print-Statistics---------
        
        //----------Stability-Statistics-----------
        fprintf(stdout, "Average Rerouting per Timestep:\n");
        for (int y = 0; y < timestep; y++) {
            if (y == 0) {
                fprintf(stdout, "%.4f,", stability[y] / TRIALS);
            } else {
                fprintf(stdout, "%.4f,", (stability[y] - stability[y-1]) / TRIALS);
            }
        }
        fprintf(stdout, "\nMax Rerouting per Node:\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", fair[y] / TRIALS);
        }
        fprintf(stdout, "\nAverage HMD Data Rate per timestep (Gbps):\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", capacity[y] / TRIALS / 1000000000);
        }
        fprintf(stdout, "\nAverage Delay Network (ms):\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", delay[y] / TRIALS);
        }
        fprintf(stdout, "\nAverage HMD Data Rate per Node over all timestep (Gbps):\n");
        for (int y = 0; y < num; y++) {
            fprintf(stdout, "%.4f,", through_all[y] / TRIALS / 1000000000);
        }
        fprintf(stdout, "\nAverage Delay Network per Node over all timestep (ms):\n");
        for (int y = 0; y < num; y++) {
            fprintf(stdout, "%.4f,", latency_all[y] / TRIALS);
        }
        
        //---------------Trials-with-Matching-Failures----------------
        fprintf(stdout, "\nMatching Failure:\n%.2f\n", count / TRIALS * 100);
        if (count > 0) {
            fprintf(stdout, "Max Unreachable count per Node in trials with Matching Failures:\n");
            for (int y = 0; y < timestep; y++) {
                fprintf(stdout, "%.4f,", reach[y] / count);
            }
            fprintf(stdout, "\nAverage HMD Data Rate per Node in trials with Matching Failures (Gbps):\n");
            for (int y = 0; y < num; y++) {
                fprintf(stdout, "%.4f,", through[y] / count / 1000000000);
            }
            fprintf(stdout, "\nAverage Delay Network per Node in trials with Matching Failures (ms):\n");
            for (int y = 0; y < num; y++) {
                fprintf(stdout, "%.4f,", latency[y] / count);
            }
        }
        
        //--------------Maximum-Unreachability-Node-----------------
        fprintf(stdout, "\nMaximum Unreachable node over all trials:\n%f\n", count2);
        if (count2 > 0) {
            fprintf(stdout, "Average HMD Data Rate per Node in trial with Max Unreachability (Gbps):\n");
            for (int y = 0; y < num; y++) {
                fprintf(stdout, "%.4f,", through_max[y] / 1000000000);
            }
            fprintf(stdout, "\nAverage Delay Network per Node in trial with Max Unreachability (ms):\n");
            for (int y = 0; y < num; y++) {
                fprintf(stdout, "%.4f,", latency_max[y]);
            }
        }
        fprintf(stdout, "\n\n");
    }
    
    //--------Print-Statistics-------
    /*
    fprintf(stdout, "Average HMD Data Rate (Gbps):\n");
    for (int y = 0; y < num-1; y++) {
        fprintf(stdout, "%.4f,", capacity[y] / TRIALS / 1000000000);
    }
    fprintf(stdout, "\nAverage Delay (ms):\n");
    for (int y = 0; y < num-1; y++) {
        fprintf(stdout, "%.4f,", delay[y] / TRIALS);
    }
    fprintf(stdout, "\n");
     */
    
    return 0;
}
