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

#define TRIALS 1000
#define GAGO 200  //max timesteps

int main(int argc, char *argv[]) {
    if (argc != 9) {
        fprintf(stderr, "./main %%ui[Number of Nodes] \
                                %%ui[X Dimension] \
                                %%ui[Y Dimension] \
                                %%ui[AP x-coordinate] \
                                %%ui[AP y-coordinate] \
                                %%ui[AP Height] \
                                %%ui[Time Step] \
                                %%ui[Num AP]\n");
        return 1;
    }
    
    unsigned int num = atoi(argv[1]);
    int length = atoi(argv[3]);
    int width = atoi(argv[2]);
    int ap_x = atoi(argv[4]);
    int ap_y = atoi(argv[5]);
    double ap_height = atof(argv[6]);
    int timestep = atoi(argv[7]);
    int num_mirrors = atoi(argv[8]);
    int region = 4;
    if (num == 0) {
        fprintf(stderr, "Number of nodes must not be 0\n");
        return 1;
    }
    time_t t;
    srand((unsigned) time(&t));
    int group_size[6] = {1, 2, 4, 5, 8, 10};
    for (int i = 0; i < sizeof(group_size); i++) {
        //ap_height = 0;
        fprintf(stdout, "group_size: %d\n", group_size[i]);
    //for (int j = 0; j <= 40; j +=2, ap_height += 0.2) {
        //fprintf(stdout, "ap_height: %.1f\n", ap_height);
    //for (num = 10; num < 110; num += 10) {
        //fprintf(stdout, "num: %u\n", num);
    //for (int region = 1; region <= 10; region++) {
        //fprintf(stdout, "region: %d\n", region);
        //double no_re_route_blockage[GAGO] = {0};
        //double num_blocked[GAGO] = {0};
        //double blockage_occurence_g[GAGO] = {0};
        //double blockage_occurence_m[GAGO] = {0};
        //double blockage_occurence_2[GAGO] = {0};
        //double greedy_blockage[GAGO] = {0};
        //double perfect_blockage[GAGO] = {0};
        //double greedy_depth2_blockage[GAGO] = {0};
        //double stability[GAGO] = {0};
        //double blockage_occurence_s[GAGO] = {0};
        //double stable[GAGO] = {0};
        //double fair[GAGO] = {0};
        //double reach[GAGO] = {0};
        double capacity[GAGO] = {0};
        double delay [GAGO] = {0};
        //struct stat height[81];
        //init_stat(height);
        double count = 0;
        for (int x = 0; x < TRIALS; x++) {
            //int fail = 0;
            //fprintf(stderr, "x:%d\n", x);
            struct graph *graph = generate_graph_group(width, length, ap_x, ap_y, ap_height, num, group_size[i]);
            fill_group(graph, group_size[i]);
            sort_group_capacity(graph);
            //sort_group_distance(graph);
            if (graph == NULL) {
                return 1;
            }
            
            //sort_height_index(graph);
            double ret = check_blockage(graph);
            //fprintf(stderr, "blocked: %f\n", ret);
            /*
            num_blocked[0] += ret;
            if (((int)(ret + 0.5)) > 0) {
                no_re_route_blockage[0]++;
            }
            if (ret > 50 ) {
                ill_posed[0]++;
            }
            */
            
            //visualize_graph(graph);
            find_parents(graph);
            find_distance(graph);
            //sort_parent_distance(graph);
            sort_parent_capacity(graph);
            //sort_parent_height(graph);
            
            //ret = greedy_matching(graph, 0);
            
            //greedy_blockage[0] += ret;
            /*
            if (((int)(ret + 0.5)) > 0) {
                //blockage_occurence_g[0]++;
            }*/
            
            
            //ret = greedy_matching_depth2(graph, ret, 0);
            /*
            greedy_depth2_blockage[0] += ret;
            if (ret > 0) {
                blockage_occurence_2[0]++;
            }*/
            
            //ret = maximal_matching(graph, 0);
            /*
            perfect_blockage[0] += ret;
            if (ret > 0) {
                blockage_occurence_m[0]++;
            }*/
            
            ret = group_matching(graph, 0);
            
            //ret = stable_matching_close(graph);
            //ret = stable_matching(graph, 0);
            //ret = stable_matching_close_height(graph);
            //ret = stable_matching_height(graph);
            
            update_capacity_delay(graph, 0);
            capacity[0] += get_capacity(graph);
            delay[0] += get_delay(graph);
            /*
            stable[0] += ret;
            if (((int)(ret + 0.5)) > 0) {
                blockage_occurence_s[0]++;
            }
            */
            /*
            if (((int)(ret + 0.5)) > 0) {
                fail = 1;
            }
            fair[0] += get_stabl(graph);
            if (fail == 1) {
                reach[0] += get_reach(graph);
            }
            */
            for (int y = 1; y < timestep; y++) {
                //update_graph(graph); //random
                update_graph_waypoint_group(width, length, graph); //random waypoint group
                //shift_index(graph);
                //sort_stability(graph);
                //sort_reachability(graph);

                ret = update_blockage(graph);
                sort_group_parent_capacity(graph);
                /*
                num_blocked[y] += ret;
                if (((int)(ret + 0.5)) > 0) {
                    no_re_route_blockage[y]++;
                }*/
                //fprintf(stderr, "num_blocked: %f\n", ret);

                //update_parents_depth2(graph);
                //update_parents(graph, y);
                //ret = update_parents_stable(graph, y);
                ret = update_parents_group(graph, y);
                /*
                reset(graph);
                find_parents(graph);
                */
                find_distance(graph);
                sort_parent_distance(graph);
                //sort_parent_capacity(graph);
                
                //ret = update_greedy(graph, y);
                //ret = update_greedy_stable(graph, y);
                /*
                //greedy_blockage[y] += ret;
                if (((int)(ret + 0.5)) > 0) {
                    blockage_occurence_g[y]++;
                }*/
                //fprintf(stderr, "remaining: %f\n\n", ret);
                
                //ret = update_depth2(graph);
                /*
                greedy_depth2_blockage[y] += ret;
                if (ret > 0) {
                    blockage_occurence_2[y]++;
                }*/
                
                //ret = update_perfect(graph, y);
                //ret = update_perfect_stable(graph, y);
                /*
                perfect_blockage[y] += ret;
                if (ret > 0) {
                    blockage_occurence_m[y]++;
                }*/

                //ret = update_stable_stable(graph, ret, y);
                //ret = update_stable_fair(graph, ret, y);
                //ret = update_group(graph, ret, y);
                //ret = update_group_stable(graph, ret, y);

                ret = update_group_fair(graph, ret, y);
                /*
                stable[y] += ret;
                if (((int)(ret + 0.5)) > 0) {
                    blockage_occurence_s[y]++;
                }*/
                
                 //Baseline
                //ret = greedy_matching(graph, y);
                //greedy_blockage[y] += ret;
                /*
                if (((int)(ret + 0.5)) > 0) {
                    blockage_occurence_g[y]++;
                }*/
                
                 /*
                ret = greedy_matching_depth2(graph, ret, y);
                greedy_depth2_blockage[0] += ret;
                if (ret > 0) {
                    blockage_occurence_2[0]++;
                }
                */
                
                 //ret = maximal_matching(graph, y);
                /*
                 perfect_blockage[y] += ret;
                 if (ret > 0) {
                 blockage_occurence_m[y]++;
                 }
                */
                
                //ret = group_matching_fair(graph, y);
                //ret = group_matching(graph , y);
                //ret = stable_matching(graph, y);
                //ret = stable_matching_fair(graph, y);
                /*
                stable[y] += ret;
                if (((int)(ret + 0.5)) > 0) {
                    blockage_occurence_s[y]++;
                }*/
                /*
                if (((int)(ret + 0.5)) > 0) {
                    fail = 1;
                }
                stability[y] += calc_stability(graph, y);
                fair[y] += get_stabl(graph);
                if (fail == 1) {
                    reach[y] += get_reach(graph);
                }
                 */
                
                update_capacity_delay(graph, y);
                capacity[y] += get_capacity(graph);
                delay[y] += get_delay(graph);
            }
            /*
            if (fail == 1) {
                count++;
            }*/
            /*
            ret = get_reach(graph);
            if (ret > count) {
                count = ret;
            }
             */
            //save_stat(graph, height);
            destroy_resources(graph);
        }
        /*
        fprintf(stdout, "height:\n");
        for (int x = 0; x < 81; x++) {
            fprintf(stdout, "%.4f,", height[x].stability / height[x].trials);
        }
        */
        fprintf(stdout, "Average HMD Data Rate (Gbps):\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", capacity[y] / TRIALS / 1000000000);
        }
        fprintf(stdout, "\nAverage Delay (ms):\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", delay[y] / TRIALS);
        }
        /*
        fprintf(stdout, "Matching Failure: %.2f\n", count / TRIALS * 100);
        fprintf(stdout, "Max Unreachability per Node:\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", reach[y] / count);
        }
        
        fprintf(stdout, "\nMax Rerouting per Node:\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", fair[y] / TRIALS);
        }
        
        fprintf(stdout, "\nAverage Rerouting per Timestep:\n");
        for (int y = 0; y < timestep; y++) {
            if (y == 0) {
                fprintf(stdout, "%.4f,", stability[y] / TRIALS);
            } else {
                fprintf(stdout, "%.4f,", (stability[y] - stability[y-1]) / TRIALS);
            }
        }*/
        
        /*
        fprintf(stdout, "\nfailure:\n");
        for (int y = 0; y < timestep; y++) {
            fprintf(stdout, "%.4f,", blockage_occurence_g[y] * 100 / TRIALS);
        }*/
        //fprintf(stdout, "%f", count);
        
        fprintf(stdout, "\n");
    }
        
    //}
    //}
    return 0;
}
