//
//  sim.c
//  
//
//  Created by Ben Ruktantichoke on 8/9/17.
//
//

#include "sim.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define SPACE 0
#define EXPO 2
#define DEV 5.8
#define BANDWIDTH 7000000000
#define PI 3.1415926535
#define LIGHT 299700000
#define GHZ 60000000000
#define LIMIT 5598720000
#define RENDER 6.1
#define BEAM 1.01
#define FOUR 2.53
#define IMAGE 62208000
#define NET 2.4

double gaussrand() {
    static double V1, V2, S;
    static int phase = 0;
    double X;
    
    if (phase == 0) {
        do {
            double U1 = (double) rand() / RAND_MAX;
            double U2 = (double) rand() / RAND_MAX;
            
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);
        
        X = V1 * sqrt(-2 * log(S)  / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
    
    phase = 1 - phase;
    
    return X;
}

double distance(struct node *node1, struct node *node2) {
    return sqrt(pow(node2->x - node1->x, 2) + pow(node2->y - node1->y, 2) + pow(node2->height - node1->height, 2));
}

void init_node(struct node *coordinate[WIDTH][LENGTH], struct node *node, int x, int y) {
    node->person = 1;
    node->x = x;
    node->y = y;
    coordinate[x][y] = node;
    node->height = ((double) (rand() % 81)) / 100 + MIN_HEIGHT;
    node->num_child = 0;
    node->num_blockers = 0;
    node->num_parent = 0;
    node->traversed = 0;
    node->idx = 0;
    memset(node->child, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->blockers, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->parent, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->distance, 0, sizeof(double) * MAX_NODE);
    memset(node->pp, 0, sizeof(struct node *) * MAX_NODE);
    node->blocked = 0;
    node->stability = 0;
    node->checked = 0;
    node->x_dest = node->x;
    node->y_dest = node->y;
    node->timer = 0;
}

void init_node_group(int width, int length, struct node *coordinate[WIDTH][LENGTH], struct node *people, int x, int y, int group_size, int idx) {
    struct node *node = &people[idx];
    node->person = 1;
    node->x = x;
    node->y = y;
    coordinate[x][y] = node;
    node->height = ((double) (rand() % 81)) / 100 + MIN_HEIGHT;
    node->num_child = 0;
    node->num_blockers = 0;
    node->num_parent = 0;
    node->traversed = 0;
    node->idx = 0;
    memset(node->child, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->blockers, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->parent, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->distance, 0, sizeof(double) * MAX_NODE);
    memset(node->distance_group, 0, sizeof(double) * MAX_GROUP);
    memset(node->pp, 0, sizeof(struct node *) * MAX_NODE);
    node->blocked = 0;
    node->stability = 0;
    node->checked = 0;
    node->x_dest = node->x;
    node->y_dest = node->y;
    node->timer = 0;
    node->marked = 0;
    node->reachability = 0;
    node->capacity = 0;
    node->delay = 0;
    
    for (int i = 1; i < group_size; i++) {
        struct node *node2 = &people[idx + i];
        node->blockers[node->num_blockers++] = node2;
        node2->person = 0;
        
        while (1) {
            int x1 = x, y1 = y, tempx = 0, tempy = 0;
            while (tempx == 0 && tempy == 0) {
                tempx = rand() % RADIUS;
                tempy = rand() % RADIUS;
            }
            x1 = (rand()%2 == 0) ? x1 + tempx : x1 - tempx;
            y1 = (rand()%2 == 0) ? y1 + tempy : y1 - tempy;

            if ((x1 >= 0) && (x1 < length) &&
                (y1 >= 0 + SPACE) && (y1 < width)) {
                if (coordinate[x1][y1] == NULL) {
                    coordinate[x1][y1] = node2;
                    node2->x = x1;
                    node2->y = y1;
                    node2->x_dest = x1 - x;
                    node2->y_dest = y1 - y;
                    break;
                }
            }
        }

        node2->height = ((double) (rand() % 81)) / 100 + MIN_HEIGHT;
        node2->num_child = 0;
        node2->num_blockers = 0;
        node2->num_parent = 0;
        node2->traversed = 0;
        node2->idx = 0;
        memset(node2->child, 0, sizeof(struct node *) * MAX_NODE);
        memset(node2->blockers, 0, sizeof(struct node *) * MAX_NODE);
        memset(node2->parent, 0, sizeof(struct node *) * MAX_NODE);
        memset(node2->distance, 0, sizeof(double) * MAX_NODE);
        memset(node2->distance_group, 0, sizeof(double) * MAX_GROUP);
        memset(node2->pp, 0, sizeof(struct node *) * MAX_NODE);
        node2->blocked = 0;
        node2->stability = 0;
        node2->checked = 0;
        node2->timer = 0;
        node2->marked = 0;
        node2->reachability = 0;
        node2->capacity = 0;
        node2->delay = 0;
    }
}

int gen_poisson(double num, double regions) {
    double lambda = regions / num;
    double seed = ((double) rand()) / RAND_MAX;
    return (int) (-log(seed) / lambda);
}

void init_node_poisson(int length, int width, struct node *coordinate[WIDTH][LENGTH], struct node *node, int region, int reg_num) {
    int len = length / region;
    int wid = (width - SPACE) / region;
    int off_x = ((reg_num - 1) / region) * len;
    int off_y = ((reg_num - 1) % region) * wid;
    node->person = 1;
    while (1) {
        int x = (rand()% len) + off_x;
        int y = (rand()% (wid)) + SPACE + off_y;
        if (x >= length || y >= width) {
            fprintf(stderr, "outside bounds\n");
            continue;
        }
        if (coordinate[x][y] == NULL) {
            node->x = x;
            node->y = y;
            coordinate[x][y] = node;
            break;
        }
    }
    node->height = ((double) (rand() % 81)) / 100 + MIN_HEIGHT;
    node->num_child = 0;
    node->num_blockers = 0;
    node->num_parent = 0;
    node->traversed = 0;
    node->idx = 0;
    memset(node->child, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->blockers, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->parent, 0, sizeof(struct node *) * MAX_NODE);
    memset(node->distance, 0, sizeof(double) * MAX_NODE);
    memset(node->distance_group, 0, sizeof(double) * MAX_GROUP);
    memset(node->pp, 0, sizeof(struct node *) * MAX_NODE);
    node->blocked = 0;
    node->stability = 0;
    node->checked = 0;
    node->x_dest = node->x;
    node->y_dest = node->y;
    node->timer = 0;
    node->marked = 0;
}

struct graph *generate_graph_unif(int width, int length, int ap_x, int ap_y, int population, double ap_height) {
    if (ap_x >= width || ap_y >= length) {
        fprintf(stderr, "AP locations is outside cooridnate\n");
        return NULL;
    }
    
    struct graph *ret = (struct graph *) malloc(sizeof(struct graph));
    //ret->coordinate = (struct node ***) malloc(sizeof(struct node *) * length * width);
    memset(ret->coordinate, 0, sizeof(struct node *) * length * width);
    
    ret->population = population;
    ret->coordinate[ap_x][ap_y] = &(ret->AP);
    ret->num_mirrors = 0;
    
    ret->AP.person = 0;
    ret->AP.x = ap_x;
    ret->AP.y = ap_y;
    ret->AP.height = ap_height;
    ret->AP.blocked = 0;
    ret->AP.num_child = 0;
    ret->AP.num_parent = 0;
    ret->AP.num_blockers = 0;
    ret->AP.num_parent = 0;
    ret->AP.traversed = 0;
    ret->AP.stability = 0;
    ret->AP.checked = 0;
    ret->AP.x_dest = ap_x;
    ret->AP.y_dest = ap_y;
    ret->AP.timer = 0;
    ret->mirrors = NULL;
    ret->AP.idx = 0;
    
    memset(ret->AP.child, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.blockers, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.parent, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.distance, 0, sizeof(double) * MAX_NODE);
    memset(ret->AP.distance_group, 0, sizeof(double) * MAX_GROUP);
    memset(ret->AP.pp, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->rr, 0, sizeof(struct node *) * MAX_NODE);
    
    ret->people = (struct node *) malloc(sizeof(struct node) * 100);
    int num = 0;
    double dlim = sqrt(ret->population);
    int lim = sqrt(ret->population);
    if (lim != ((int) (dlim + 0.5))) {
        lim += 1;
    }
    int intx = width/sqrt(ret->population);
    int inty = length/sqrt(ret->population);
    for (int i = intx; i < width; i+= intx) {
        for (int j = inty; j < length; j+= inty) {
            init_node(ret->coordinate, &ret->people[num++], i, j);
            if (num >= ret->population) {
                return ret;
            }
        }
    }
    
    return ret;
}

struct graph *generate_graph_group(int width, int length, int ap_x, int ap_y, double ap_height, int population, int group_size) {
    if (ap_x >= width || ap_y >= length) {
        fprintf(stderr, "AP location is outside cooridnate\n");
        return NULL;
    }
    
    struct graph *ret = (struct graph *) malloc(sizeof(struct graph));
    //ret->coordinate = (struct node ***) malloc(sizeof(struct node *) * length * width);
    memset(ret->coordinate, 0, sizeof(struct node *) * length * width);
    
    ret->population = population;
    ret->coordinate[ap_x][ap_y] = &(ret->AP);
    ret->num_mirrors = 0;
    
    ret->AP.person = 0;
    ret->AP.x = ap_x;
    ret->AP.y = ap_y;
    ret->AP.height = ap_height;
    ret->AP.blocked = 0;
    ret->AP.num_child = 0;
    ret->AP.num_parent = 0;
    ret->AP.num_blockers = 0;
    ret->AP.num_parent = 0;
    ret->AP.traversed = 0;
    ret->AP.stability = 0;
    ret->AP.checked = 0;
    ret->AP.x_dest = ap_x;
    ret->AP.y_dest = ap_y;
    ret->AP.timer = 0;
    ret->mirrors = NULL;
    ret->AP.idx = 0;
    ret->AP.marked = 0;
    ret->AP.reachability = 0;
    ret->AP.capacity = 0;
    ret->AP.delay = 0;
    
    memset(ret->AP.child, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.blockers, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.parent, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.distance, 0, sizeof(double) * MAX_NODE);
    memset(ret->AP.distance_group, 0, sizeof(double) * MAX_GROUP);
    memset(ret->AP.pp, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->rr, 0 , sizeof(struct node *) * MAX_NODE);
    
    ret->people = (struct node *) malloc(sizeof(struct node) * population);
    int num = 0;
    int i = 0;
    int l = ret->population / group_size;
    for (; i < l; i++) {
        while (1) {
            int x = (rand() % (length - 10)) + 5;
            int y = (rand() % (width - SPACE - 10)) + 5;
            if (ret->coordinate[x][y + SPACE] == NULL) {
                init_node_group(width, length, ret->coordinate, ret->people, x, y + SPACE, group_size, num);
                num += group_size;
                break;
            }
        }
    }
    //fprintf(stderr, "leftover: %d\n",ret->population-(i)*group_size);
    if (i * group_size < ret->population) {
        while (1) {
            int x = (rand() % (length - 10)) + 5;
            int y = (rand() % (width - SPACE - 10)) + 5;
            if (ret->coordinate[x][y + SPACE] == NULL) {
                init_node_group(width, length, ret->coordinate, ret->people, x, y + SPACE, ret->population-i*group_size, num);
                break;
            }
        }
    }
    
    for (int i = 0; i < ret->population; i++) {
        ret->rr[i] = &ret->people[i];
    }
 
    return ret;
}

void fill_group(struct graph *graph, int group_size) {
    if (group_size > 1) {
        int i = 0;
        for (; i + group_size < graph->population; i += group_size) {
            for (int j = i; j < i + group_size; j++) {
                if (graph->people[j].person == 1) {
                    continue;
                }
                for (int w = i; w < i + group_size; w++) {
                    if (w == j) {
                        continue;
                    }
                    graph->people[j].blockers[graph->people[j].num_blockers++] = &graph->people[w];
                }
            }
        }
        if (i > graph->population) {
            for (i -= group_size; i < graph->population; i++) {
                if (graph->people[i].person == 1) {
                    continue;
                }
                for (int w = i; w < graph->population; w++) {
                    if (w == i) {
                        continue;
                    }
                    graph->people[i].blockers[graph->people[i].num_blockers++] = &graph->people[w];
                }
            }
        }
    }
    
    for (int i = 0; i < graph->population; i++) {
        for (int j = 0; j < graph->people[i].num_blockers; j++) {
            graph->people[i].distance_group[j] = distance(&graph->people[i], graph->people[i].blockers[j]);
        }
    }
}

void sort_group_distance(struct graph *graph) {
    for (int i = 0; i < graph->population; i++) {
        struct node *temp_node = NULL;
        double temp;
        for (int y = 0; y < graph->people[i].num_blockers; y++) {
            for (int z = y + 1; z < graph->people[i].num_blockers; z++) {
                if (graph->people[i].distance_group[y] > graph->people[i].distance_group[z]) {
                    temp = graph->people[i].distance_group[y];
                    graph->people[i].distance_group[y] = graph->people[i].distance_group[z];
                    graph->people[i].distance_group[z] = temp;
                    temp_node = graph->people[i].blockers[y];
                    graph->people[i].blockers[y] = graph->people[i].blockers[z];
                    graph->people[i].blockers[z] = temp_node;
                } else if (graph->people[i].distance_group[y] == graph->people[i].distance_group[z]) {
                    if (graph->people[i].blockers[y]->height < graph->people[i].blockers[z]->height) {
                        temp = graph->people[i].distance_group[y];
                        graph->people[i].distance_group[y] = graph->people[i].distance_group[z];
                        graph->people[i].distance_group[z] = temp;
                        temp_node = graph->people[i].blockers[y];
                        graph->people[i].blockers[y] = graph->people[i].blockers[z];
                        graph->people[i].blockers[z] = temp_node;
                    }
                }
            }
        }
    }
}

void sort_group_capacity(struct graph *graph) {
    for (int i = 0; i < graph->population; i++) {
        struct node *temp_node = NULL;
        double temp;
        for (int y = 0; y < graph->people[i].num_blockers; y++) {
            for (int z = y + 1; z < graph->people[i].num_blockers; z++) {
                if (graph->people[i].blockers[y].capacity > graph->people[i].blockers[z].capacity) {
                    temp = graph->people[i].distance_group[y];
                    graph->people[i].distance_group[y] = graph->people[i].distance_group[z];
                    graph->people[i].distance_group[z] = temp;
                    temp_node = graph->people[i].blockers[y];
                    graph->people[i].blockers[y] = graph->people[i].blockers[z];
                    graph->people[i].blockers[z] = temp_node;
                } else if (graph->people[i].blockers[y].capacity == graph->people[i].blockers[z].capacity) {
                    if (graph->people[i].distance_group[y] < graph->people[i].distance_group[z]) {
                        temp = graph->people[i].distance_group[y];
                        graph->people[i].distance_group[y] = graph->people[i].distance_group[z];
                        graph->people[i].distance_group[z] = temp;
                        temp_node = graph->people[i].blockers[y];
                        graph->people[i].blockers[y] = graph->people[i].blockers[z];
                        graph->people[i].blockers[z] = temp_node;
                    }
                }
            }
        }
    }
}

struct graph *generate_graph_poisson(int width, int length, int ap_x, int ap_y, int num, double ap_height, int region) {
    if (num > length * width) {
        fprintf(stderr, "Number of nodes requested is larger than coordinate space\n");
        return NULL;
    }
    if (ap_x >= width || ap_y >= length) {
        fprintf(stderr, "AP locations is outside cooridnate\n");
        return NULL;
    }
    
    struct graph *ret = (struct graph *) malloc(sizeof(struct graph));
    //ret->coordinate = (struct node ***) malloc(sizeof(struct node *) * length * width);
    memset(ret->coordinate, 0, sizeof(struct node *) * length * width);
    
    ret->population = num;
    ret->coordinate[ap_x][ap_y] = &(ret->AP);
    ret->num_mirrors = 0;
    
    ret->AP.person = 0;
    ret->AP.x = ap_x;
    ret->AP.y = ap_y;
    ret->AP.height = ap_height;
    ret->AP.blocked = 0;
    ret->AP.num_child = 0;
    ret->AP.num_parent = 0;
    ret->AP.num_blockers = 0;
    ret->AP.num_parent = 0;
    ret->AP.traversed = 0;
    ret->AP.stability = 0;
    ret->AP.checked = 0;
    ret->AP.x_dest = ap_x;
    ret->AP.y_dest = ap_y;
    ret->AP.timer = 0;
    ret->AP.idx = 0;
    ret->AP.marked = 0;
    ret->mirrors = NULL;
    
    memset(ret->AP.child, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.blockers, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.parent, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->AP.distance, 0, sizeof(double) * MAX_NODE);
    memset(ret->AP.distance_group, 0, sizeof(double) * MAX_GROUP);
    memset(ret->AP.pp, 0, sizeof(struct node *) * MAX_NODE);
    memset(ret->rr, 0, sizeof(struct node *) * MAX_NODE);
    
    ret->people = (struct node *) malloc(sizeof(struct node) * num);
    
    int regions = region * region;
    int accum = num;
    int node_id = 0;
    while (accum > 0) {
        for (int i = 1; i <= regions; i++) {
            int poisson = gen_poisson((double)num, (double)regions);
            if (accum - poisson > 0) {
                accum -= poisson;
            } else {
                poisson = accum;
                accum -= poisson;
            }
            for (int j = 0; j < poisson; j++) {
                init_node_poisson(width, length, ret->coordinate, &ret->people[node_id++], region, i);
            }
        }
    }
    
    return ret;
}

void visualize_graph(int width, int length, struct graph *graph) {
    for (int j = length - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            if (graph->coordinate[i][j] != NULL) {
                //fprintf(stderr, "blocked:%d\n", graph->coordinate[i][j]->blocked);
                if (graph->coordinate[i][j] == &(graph->AP)) {
                    fprintf(stderr, "* ");
                } /*else if (graph->coordinate[i][j]->person == 0) {
                    fprintf(stderr, "# ");
                }*/ else if (graph->coordinate[i][j]->blocked > 0) {
                    fprintf(stderr, "+ ");
                } else {
                    fprintf(stderr, "O ");
                }
            } else {
                fprintf(stderr, "- ");
            }
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
    usleep(1000000);
}

void visualize_stability(int width, int length, struct graph *graph) {
    for (int j = length - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            if (graph->coordinate[i][j] != NULL) {
                //fprintf(stderr, "blocked:%d\n", graph->coordinate[i][j]->blocked);
                if (graph->coordinate[i][j] == &(graph->AP)) {
                    fprintf(stderr, "*  ");
                } /*else if (graph->coordinate[i][j]->person == 0) {
                   fprintf(stderr, "# ");
                   }*/ else {
                       if ((int) (graph->coordinate[i][j]->stability + 0.5) > 9) {
                           fprintf(stderr, "%d ",(int) (graph->coordinate[i][j]->stability + 0.5));
                       } else {
                           fprintf(stderr, "%d  ",(int) (graph->coordinate[i][j]->stability + 0.5));
                       }
                   }
            } else {
                fprintf(stderr, "-  ");
            }
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
    usleep(1000000);
}

void visualize_reachability(int width, int length, struct graph *graph) {
    for (int j = length - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            if (graph->coordinate[i][j] != NULL) {
                //fprintf(stderr, "blocked:%d\n", graph->coordinate[i][j]->blocked);
                if (graph->coordinate[i][j] == &(graph->AP)) {
                    fprintf(stderr, "* ");
                } /*else if (graph->coordinate[i][j]->person == 0) {
                   fprintf(stderr, "# ");
                   }*/ else if (graph->coordinate[i][j]->blocked == 1) {
                       fprintf(stderr, "+ ");
                   } else {
                       fprintf(stderr, "O ");
                   }
            } else {
                fprintf(stderr, "- ");
            }
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
    usleep(1000000);
}

void destroy_resources(struct graph *graph) {
    //free(graph->coordinate);
    free(graph->people);
    if (graph->num_mirrors != 0) {
        free(graph->mirrors);
    }
    free(graph);
}

double DtoLOS(struct node *AP, struct node *target, struct node *candidate) {
    double ret;
    double num;
    double denom;
    num = fabs((target->y - AP->y)*candidate->x -
               (target->x - AP->x)*candidate->y +
               target->x*AP->y -
               target->y*AP->x);
    denom = sqrt(pow(target->y - AP->y, 2) + pow(target->x - AP->x, 2));
    ret = num/denom;
    return ret;
}

double HatLOS(struct node *AP, struct node *target, struct node *candidate) {
    double ret;
    double u_x = target->x - AP->x;
    double u_y = target->y - AP->y;
    double u_z = target->height - AP->height;
    double pq_x = candidate->x - AP->x;
    double pq_y = candidate->y - AP->y;
    double pq_z = candidate->height - AP->height;
    
    double dot = u_x * pq_x + u_y * pq_y + u_z * pq_z;
    double temp = dot / (pow(u_x, 2) +
                         pow(u_y, 2) +
                         pow(u_z, 2));

    double du_z = u_z * temp;
    double w_z = pq_z - du_z;
    
    ret = candidate->height - w_z;
    return ret;
}

//return 1 if blocked, 0 if LOS exists
int check_blockage_node(struct node *node, struct node *target, struct node *coordinate[WIDTH][LENGTH], int chance) {
    int temp_nx = (int) (node->x + 0.5);
    int temp_tx = (int) (target->x + 0.5);
    int temp_ny = (int) (node->y + 0.5);
    int temp_ty = (int) (target->y + 0.5);
    //fprintf(stderr, "temp_nx: %d, temp_tx %d\ntemp_ny: %d, temp_ty: %d\n", temp_nx, temp_tx, temp_ny, temp_ty);
    if (target->x > node->x) {
        for (int x = temp_nx; x <= temp_tx; x++) {
            if (target->y > node->y) {
                for (int y = temp_ny; y <= temp_tx; y++) {
                    //fprintf(stderr, "x: %d, y: %d\n", x, y);
                    if (coordinate[x][y] != NULL) {
                        if ((x == temp_tx && y == temp_ty) ||
                            (x == temp_nx && y == temp_ny)) {
                            continue;
                        }
                        
                        double temp = DtoLOS(target, node, coordinate[x][y]);
                        if (temp < FAT) {
                            if (coordinate[x][y]->height >= target->height &&
                                coordinate[x][y]->height >= node->height) { //Candidate Taller than AP and Target
                                return 1;
                            } else if (coordinate[x][y]->height >= target->height ||
                                       coordinate[x][y]->height >= node->height) { //Somewhere in between
                                if (HatLOS(target, node, coordinate[x][y]) <
                                    coordinate[x][y]->height) {
                                    return 1;
                                }
                            } else { //Candidate shorter than AP and Target
                                continue;
                            }
                        }
                    }
                }
            } else {
                for (int y = temp_ty; y <= temp_ny; y++) {
                    //fprintf(stderr, "x: %d, y: %d\n", x, y);
                    if (coordinate[x][y] != NULL) {
                        if ((x == temp_tx && y == temp_ty) ||
                            (x == temp_nx && y == temp_ny)) {
                            continue;
                        }

                        double temp = DtoLOS(target, node, coordinate[x][y]);
                        if (temp < FAT) {
                            if (coordinate[x][y]->height >= target->height &&
                                coordinate[x][y]->height >= node->height) {
                                return 1;
                            } else if (coordinate[x][y]->height >= target->height ||
                                       coordinate[x][y]->height >= node->height) {
                                if (HatLOS(target, node, coordinate[x][y]) <
                                    coordinate[x][y]->height) {
                                    return 1;
                                }
                            } else {
                                continue;
                            }
                        }
                    }
                }
            }
        }
    } else {
        for (int x = temp_tx; x <= temp_nx; x++) {
            if (target->y > node->y) {
                for (int y = temp_ny; y <= temp_ty; y++) {
                    //fprintf(stderr, "x: %d, y: %d\n", x, y);
                    if (coordinate[x][y] != NULL) {
                        if ((x == temp_tx && y == temp_ty) ||
                            (x == temp_nx && y == temp_ny)) {
                            continue;
                        }
                        
                        double temp = DtoLOS(target, node, coordinate[x][y]);
                        if (temp < FAT) {
                            if (coordinate[x][y]->height >= target->height &&
                                coordinate[x][y]->height >= node->height) {
                                return 1;
                            } else if (coordinate[x][y]->height >= target->height ||
                                       coordinate[x][y]->height >= node->height) {
                                if (HatLOS(target, node, coordinate[x][y]) <
                                    coordinate[x][y]->height) {
                                    return 1;
                                }
                            } else {
                                continue;
                            }
                        }
                    }
                }
            } else {
                for (int y = temp_ty; y <= temp_ny; y++) {
                    //fprintf(stderr, "x: %d, y: %d\n", x, y);
                    if (coordinate[x][y] != NULL) {
                        if ((x == temp_tx && y == temp_ty) ||
                            (x == temp_nx && y == temp_ny)) {
                            continue;
                        }
                        
                        double temp = DtoLOS(target, node, coordinate[x][y]);
                        if (temp < FAT) {
                            if (coordinate[x][y]->height >= target->height &&
                                coordinate[x][y]->height >= node->height) {
                                return 1;
                            } else if (coordinate[x][y]->height >= target->height ||
                                       coordinate[x][y]->height >= node->height) {
                                if (HatLOS(target, node, coordinate[x][y]) <
                                    coordinate[x][y]->height) {
                                    return 1;
                                }
                            } else {
                                continue;
                            }
                        }
                    }
                }
            }
        }
    }
    
    int l;
    //10% Other obstacles block chance
    switch (chance) {
        case 0:
            return 0;
        case 1:
            l = rand() % 9;
            if (l == 0) {
                return 1;
            }
        default:
            return 0;
    }
}

double calc_capacity(struct node *node, struct node *target) {
    double d = distance(node, target);
    return BANDWIDTH * log2(1+pow(10, (116 - EXPO * 10 * log10(4 * PI * d * GHZ / LIGHT) + DEV * gaussrand())/10));
}

//Construct list of blocked nodes in graph->AP.blockers and LOS nodes in graph->AP.child
double check_blockage(struct graph *graph) {
    for (int x = 0; x < graph->population; x++) {
        if (check_blockage_node(graph->rr[x], &graph->AP, graph->coordinate, 0) == 1) {
            graph->rr[x]->blocked = 1;
            graph->rr[x]->capacity = 0;
            graph->AP.blockers[graph->AP.num_blockers++] = graph->rr[x];
        } else {
            graph->rr[x]->capacity = calc_capacity(graph->rr[x], &graph->AP);
            //fprintf(stderr, "capacity: %f\n", graph->rr[x]->capacity);
            graph->AP.child[graph->AP.num_child++] = graph->rr[x];
        }
    }

    return (double) graph->AP.num_blockers;
}

//used only for Depth2
double update_blockage_d2(struct graph *graph) {
    memset(graph->AP.blockers, 0, sizeof(struct node *) * graph->AP.num_blockers);
    graph->AP.num_blockers = 0;
    memset(graph->AP.child, 0, sizeof(struct node *) * graph->AP.num_child);
    graph->AP.num_child = 0;
    
    for (int x = 0; x < graph->population; x++) {
        if (check_blockage_node(graph->rr[x], &graph->AP, graph->coordinate, 0) == 1) {
            if (graph->rr[x]->blocked == 0 &&
                graph->rr[x]->num_parent == 0) {
                graph->rr[x]->stability++;
            }
            graph->rr[x]->blocked = 1;
            graph->AP.blockers[graph->AP.num_blockers++] = graph->rr[x];
        } else {
            graph->rr[x]->blocked = 0;
            graph->AP.child[graph->AP.num_child++] = graph->rr[x];
        }
    }
    
    return graph->AP.num_blockers;
}

double update_blockage(struct graph *graph) {
    memset(graph->AP.blockers, 0, sizeof(struct node *) * graph->AP.num_blockers);
    graph->AP.num_blockers = 0;
    memset(graph->AP.child, 0, sizeof(struct node *) * graph->AP.num_child);
    graph->AP.num_child = 0;
    
    for (int x = 0; x < graph->population; x++) {
        if (check_blockage_node(graph->rr[x], &graph->AP, graph->coordinate, 0) == 1) {
            graph->rr[x]->blocked = 1;
            graph->rr[x]->capacity = 0;
            graph->rr[x]->delay = 0;
            graph->AP.blockers[graph->AP.num_blockers++] = graph->rr[x];
        } else {
            graph->rr[x]->blocked = 0;
            graph->rr[x]->capacity = calc_capacity(graph->rr[x], &graph->AP);
            //fprintf(stderr, "capacity: %f\n", graph->rr[x]->capacity);
            graph->AP.child[graph->AP.num_child++] = graph->rr[x];
        }
    }
    
    return graph->AP.num_blockers;
}

//Find potential parents for each blocked node
void find_parents(struct graph *graph) {
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        for (int y = 0; y < graph->AP.num_child; y++) {
            if (check_blockage_node(graph->AP.blockers[x], graph->AP.child[y], graph->coordinate, 0) == 0) {
                graph->AP.blockers[x]->parent[graph->AP.blockers[x]->num_parent++] = graph->AP.child[y];
            }
        }
    }
}

void find_distance(struct graph *graph) {
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
            graph->AP.blockers[x]->distance[y] = distance(graph->AP.blockers[x], graph->AP.blockers[x]->parent[y]);
        }
    }
}

void sort_parent_distance(struct graph *graph) {
    struct node *temp_node = NULL;
    double temp;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
            for (int z = y + 1; z < graph->AP.blockers[x]->num_parent; z++) {
                if (graph->AP.blockers[x]->distance[y] > graph->AP.blockers[x]->distance[z]) {
                    temp = graph->AP.blockers[x]->distance[y];
                    graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                    graph->AP.blockers[x]->distance[z] = temp;
                    temp_node = graph->AP.blockers[x]->parent[y];
                    graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                    graph->AP.blockers[x]->parent[z] = temp_node;
                } else if (graph->AP.blockers[x]->distance[y] == graph->AP.blockers[x]->distance[z]) {
                    if (graph->AP.blockers[x]->parent[y]->height < graph->AP.blockers[x]->parent[z]->height) {
                        temp_node = graph->AP.blockers[x]->parent[y];
                        graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                        graph->AP.blockers[x]->parent[z] = temp_node;
                        temp = graph->AP.blockers[x]->distance[y];
                        graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                        graph->AP.blockers[x]->distance[z] = temp;
                    }
                }
            }
        }
    }
}

void sort_parent_capacity(struct graph *graph) {
    struct node *temp_node = NULL;
    double temp;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
            for (int z = y + 1; z < graph->AP.blockers[x]->num_parent; z++) {
                if (graph->AP.blockers[x]->parent[y]->capacity > graph->AP.blockers[x]->parent[z]->capacity) {
                    temp = graph->AP.blockers[x]->distance[y];
                    graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                    graph->AP.blockers[x]->distance[z] = temp;
                    temp_node = graph->AP.blockers[x]->parent[y];
                    graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                    graph->AP.blockers[x]->parent[z] = temp_node;
                } else if (graph->AP.blockers[x]->parent[y]->capacity == graph->AP.blockers[x]->parent[z]->capacity) {
                    if (graph->AP.blockers[x]->distance[y] < graph->AP.blockers[x]->distance[z]) {
                        temp_node = graph->AP.blockers[x]->parent[y];
                        graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                        graph->AP.blockers[x]->parent[z] = temp_node;
                        temp = graph->AP.blockers[x]->distance[y];
                        graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                        graph->AP.blockers[x]->distance[z] = temp;
                    }
                }
            }
        }
    }
}

void sort_group_parent_capacity(struct graph *graph) {
    for (int i = 0; i < graph->population; i++) {
        struct node *temp_node = NULL;
        double temp;
        for (int y = 0; y < graph->people[i].num_blockers; y++) {
            for (int z = y + 1; z < graph->people[i].num_blockers; z++) {
                if (graph->people[i].blockers[y].capacity > graph->people[i].blockers[z].capacity) {
                    if (graph->people[i].traversed == 1) {
                        if (graph->people[x].idx == y) {
                            graph->people[x].idx = z;
                        } else if (graph->people[x].idx == z) {
                            graph->people[x].idx = y;
                        }
                    }
                    temp = graph->people[i].distance_group[y];
                    graph->people[i].distance_group[y] = graph->people[i].distance_group[z];
                    graph->people[i].distance_group[z] = temp;
                    temp_node = graph->people[i].blockers[y];
                    graph->people[i].blockers[y] = graph->people[i].blockers[z];
                    graph->people[i].blockers[z] = temp_node;
                } else if (graph->people[i].blockers[y].capacity == graph->people[i].blockers[z].capacity) {
                    if (graph->people[i].distance_group[y] < graph->people[i].distance_group[z]) {
                        if (graph->people[i].traversed == 1) {
                            if (graph->people[x].idx == y) {
                                graph->people[x].idx = z;
                            } else if (graph->people[x].idx == z) {
                                graph->people[x].idx = y;
                            }
                        }
                        temp = graph->people[i].distance_group[y];
                        graph->people[i].distance_group[y] = graph->people[i].distance_group[z];
                        graph->people[i].distance_group[z] = temp;
                        temp_node = graph->people[i].blockers[y];
                        graph->people[i].blockers[y] = graph->people[i].blockers[z];
                        graph->people[i].blockers[z] = temp_node;
                    }
                }
            }
        }
    }
}

void sort_parent_height(struct graph *graph) {
    double temp;
    struct node *temp_node = NULL;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
            for (int z = y + 1; z < graph->AP.blockers[x]->num_parent; z++) {
                if (graph->AP.blockers[x]->parent[y]->height < graph->AP.blockers[x]->parent[z]->height) {
                    temp_node = graph->AP.blockers[x]->parent[y];
                    graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                    graph->AP.blockers[x]->parent[z] = temp_node;
                    temp = graph->AP.blockers[x]->distance[y];
                    graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                    graph->AP.blockers[x]->distance[z] = temp;
                } else if (graph->AP.blockers[x]->parent[y]->height == graph->AP.blockers[x]->parent[z]->height) {
                    if (graph->AP.blockers[x]->distance[y] > graph->AP.blockers[x]->distance[z]) {
                        temp = graph->AP.blockers[x]->distance[y];
                        graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                        graph->AP.blockers[x]->distance[z] = temp;
                        temp_node = graph->AP.blockers[x]->parent[y];
                        graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                        graph->AP.blockers[x]->parent[z] = temp_node;
                    }
                }
            }
        }
    }
}

void update_parents(struct graph *graph, int t) {
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        memset(graph->AP.blockers[x]->parent, 0, sizeof(struct node *) * MAX_NODE);
        memset(graph->AP.blockers[x]->distance, 0, sizeof(double) * MAX_NODE);
        graph->AP.blockers[x]->num_parent = 0;
        for (int y = 0; y < graph->AP.num_child; y++) {
            if (check_blockage_node(graph->AP.blockers[x], graph->AP.child[y], graph->coordinate, 0) == 0) {
                graph->AP.blockers[x]->parent[graph->AP.blockers[x]->num_parent++] = graph->AP.child[y];
            }
        }
        
        if (graph->AP.blockers[x]->num_child == 1) {
            graph->AP.blockers[x]->num_child = 0;
            graph->AP.blockers[x]->child[0] = NULL;
        }
    }
    
    for (int x = 0; x < graph->AP.num_child; x++) {
        if (graph->AP.child[x]->num_child == 1) {
            if (check_blockage_node(graph->AP.child[x],graph->AP.child[x]->child[0], graph->coordinate, 0) == 1) {
                graph->AP.child[x]->num_child = 0;
                graph->AP.child[x]->child[0] = NULL;
            } else {
                graph->AP.child[x]->child[0]->pp[t] = graph->AP.child[x];
                graph->AP.child[x]->child[0]->checked = 1;
            }
        }
    }
}

void update_parents_depth2(struct graph *graph) {
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->num_parent == 0) { //If blocked node did not previously have a parent, then it must previously be a LOS node
            if (graph->AP.blockers[x]->num_child == 1) {//LOS nodes must have their childs reset
                //fprintf(stderr, "child: %p\n", graph->AP.blockers[x]->child[0]);
                if (graph->AP.blockers[x]->child[0]->num_child == 1) {//handling depth 2 child
                    //fprintf(stderr, "child->child: %p\n", graph->AP.blockers[x]->child[0]->child[0]);
                    graph->AP.blockers[x]->child[0]->num_child = 0;
                    graph->AP.blockers[x]->child[0]->child[0]->stability++;
                    if (graph->AP.blockers[x]->child[0]->child[0]->blocked == 0) {//If node isn't blocked then it switches back to LOS
                        graph->AP.blockers[x]->child[0]->child[0]->num_parent = 0;
                        memset(graph->AP.blockers[x]->child[0]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                    }
                    graph->AP.blockers[x]->child[0]->child[0] = NULL;
                }
                //fprintf(stderr, "mid\n");
                graph->AP.blockers[x]->num_child = 0;
                graph->AP.blockers[x]->child[0]->stability++;
                if (graph->AP.blockers[x]->child[0]->blocked == 0) {//Switch back to LOS
                    graph->AP.blockers[x]->child[0]->num_parent = 0;
                    memset(graph->AP.blockers[x]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                }
                graph->AP.blockers[x]->child[0] = NULL;
                //fprintf(stderr, "end\n");
            }
        }
        //Reset possible parent list
        memset(graph->AP.blockers[x]->parent, 0, sizeof(struct node *) * MAX_NODE);
        graph->AP.blockers[x]->num_parent = 0;
        for (int y = 0; y < graph->AP.num_child; y++) {
            if (check_blockage_node(graph->AP.blockers[x], graph->AP.child[y], graph->coordinate, 0) == 0) {
                graph->AP.blockers[x]->parent[graph->AP.blockers[x]->num_parent++] = graph->AP.child[y];
            }
        }
    }
    
    //Now we want to find all nodes that should use same link. This happens only if the same links have LOS to the AP
    for (int x = 0; x < graph->population; x++) {
        if (graph->people[x].marked == 0) {
            graph->people[x].marked = 1;
            if (graph->people[x].num_child == 1) {
                if (check_blockage_node(&graph->people[x], graph->people[x].child[0], graph->coordinate, 0) == 1) {
                    if (graph->people[x].child[0]->num_child == 1) {
                        graph->people[x].child[0]->child[0]->stability++;
                        graph->people[x].child[0]->num_child = 0;
                        graph->people[x].child[0]->child[0]->checked = 0;
                        graph->people[x].child[0]->child[0]->traversed = 0;
                        if (graph->people[x].child[0]->child[0]->blocked == 0) {
                            graph->people[x].child[0]->child[0]->num_parent = 0;
                            memset(graph->people[x].child[0]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                        }
                        graph->people[x].child[0]->child[0]->marked = 1;
                        graph->people[x].child[0]->child[0] = NULL;
                    }
                    graph->people[x].child[0]->stability++;
                    graph->people[x].child[0]->checked = 0;
                    graph->people[x].child[0]->traversed = 0;
                    if (graph->people[x].child[0]->blocked == 0) {
                        graph->people[x].child[0]->num_parent = 0;
                        memset(graph->people[x].child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                    }
                    graph->people[x].num_child = 0;
                    graph->people[x].marked = 1;
                    graph->people[x].child[0] = NULL;
                } else {
                    graph->people[x].child[0]->checked = 1;
                    graph->people[x].child[0]->traversed = 1;
                    graph->people[x].child[0]->marked = 1;
                    if (graph->people[x].child[0]->num_child == 1) {
                        if (graph->people[x].blocked == 1 || check_blockage_node(graph->people[x].child[0], graph->people[x].child[0]->child[0], graph->coordinate, 0) == 1) {
                            graph->people[x].child[0]->child[0]->stability++;
                            graph->people[x].child[0]->num_child = 0;
                            graph->people[x].child[0]->child[0]->checked = 0;
                            graph->people[x].child[0]->child[0]->traversed = 0;
                            if (graph->people[x].child[0]->child[0]->blocked == 0) {
                                graph->people[x].child[0]->child[0]->num_parent = 0;
                                memset(graph->people[x].child[0]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                            }
                            graph->people[x].child[0]->child[0]->marked = 1;
                            graph->people[x].child[0]->child[0] = NULL;
                        } else {
                            graph->people[x].child[0]->child[0]->marked = 1;
                            graph->people[x].child[0]->child[0]->checked = 1;
                            graph->people[x].child[0]->child[0]->traversed = 2;
                        }
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->marked = 0;
        if (graph->AP.blockers[x]->checked == 0) {//quick fix
            if (graph->AP.blockers[x]->num_child == 1) {
                if (graph->AP.blockers[x]->child[0]->blocked == 0) {
                    if (graph->AP.blockers[x]->child[0]->num_parent > 0) {
                        graph->AP.blockers[x]->child[0]->num_parent = 0;
                        memset(graph->AP.blockers[x]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                    }
                }
                graph->AP.blockers[x]->child[0]->checked = 0;
                graph->AP.blockers[x]->child[0]->traversed = 0;
                graph->AP.blockers[x]->num_child = 0;
                graph->AP.blockers[x]->child[0] = NULL;
            }
        }
    }
}

double update_parents_group(struct graph *graph, int t) {
    double ret = 0;
    for (int x = 0; x < graph->AP.num_child; x++) {
        if (graph->AP.child[x]->num_child > 0) {
            if (check_blockage_node(graph->AP.child[x],
                                    graph->AP.child[x]->child[0],
                                    graph->coordinate, 0) == 1) {
                //Modify blocked child
                graph->AP.child[x]->child[0]->idx = 0;
                graph->AP.child[x]->child[0]->traversed = 0;
                
                //Modify Node
                graph->AP.child[x]->num_child = 0;
                graph->AP.child[x]->child[0] = NULL;
            } else {
                graph->AP.child[x]->child[0]->marked = 1;
                graph->AP.child[x]->child[0]->pp[t] = graph->AP.child[x];
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->num_child > 0) {//Modify blocked child
            graph->AP.blockers[x]->child[0]->idx = 0;
            graph->AP.blockers[x]->child[0]->traversed = 0;
            graph->AP.blockers[x]->num_child = 0;
            graph->AP.blockers[x]->child[0] = NULL;
        }
        
        if (graph->AP.blockers[x]->marked == 0) {
            graph->AP.blockers[x]->idx = 0;
            graph->AP.blockers[x]->num_parent = 0;
            memset(graph->AP.blockers[x]->parent, 0, sizeof(struct node *) * MAX_NODE);
            memset(graph->AP.blockers[x]->distance, 0, sizeof(double) * MAX_NODE);
            
            for (int y = 0; y < graph->AP.num_child; y++) { //Find parents
                if (check_blockage_node(graph->AP.blockers[x], graph->AP.child[y], graph->coordinate, 0) == 0) {
                    graph->AP.blockers[x]->parent[graph->AP.blockers[x]->num_parent++] = graph->AP.child[y];
                }
            }
            
        } else {
            ret++;
            graph->AP.blockers[x]->blocked = 2;
        }
        
    }
    
    return ret;
}

//returns number of blocked nodes that can use same link
double update_parents_stable(struct graph *graph, int t) {
    double ret = 0;
    for (int x = 0; x < graph->AP.num_child; x++) {
        if (graph->AP.child[x]->num_child == 1) {
            if (check_blockage_node(graph->AP.child[x],
                                    graph->AP.child[x]->child[0],
                                    graph->coordinate, 0) == 1) {
                graph->AP.child[x]->num_child = 0;
                graph->AP.child[x]->child[0] = NULL;
            } else {
                if (graph->AP.child[x]->child[0]->blocked == 1) {
                    ret++;
                    graph->AP.child[x]->child[0]->blocked = 2;
                }
                graph->AP.child[x]->child[0]->checked = 1;
                graph->AP.child[x]->child[0]->pp[t] = graph->AP.child[x];
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->num_child > 0) {
            graph->AP.blockers[x]->num_child = 0;
            graph->AP.blockers[x]->child[0] = NULL;
        }
        graph->AP.blockers[x]->idx = 0;
        graph->AP.blockers[x]->num_parent = 0;
        memset(graph->AP.blockers[x]->parent, 0, sizeof(struct node *) * MAX_NODE);
        memset(graph->AP.blockers[x]->distance, 0, sizeof(double) * MAX_NODE);
        
        for (int y = 0; y < graph->AP.num_child; y++) {
            if (check_blockage_node(graph->AP.blockers[x], graph->AP.child[y], graph->coordinate, 0) == 0) {
                graph->AP.blockers[x]->parent[graph->AP.blockers[x]->num_parent++] = graph->AP.child[y];
            }
        }
    }

    return ret;
}

double update_parents_stable_height(struct graph *graph, int t) {
    double ret = 0;
    for (int x = 0; x < graph->AP.num_child; x++) {
        if (graph->AP.child[x]->num_child == 1) {
            if (check_blockage_node(graph->AP.child[x],
                                    graph->AP.child[x]->child[0],
                                    graph->coordinate, 0) == 1) {
                graph->AP.child[x]->num_child = 0;
                graph->AP.child[x]->child[0] = NULL;
            } else {
                if (graph->AP.child[x]->child[0]->blocked == 1) {
                    ret++;
                    graph->AP.blockers[x]->blocked = 2;
                }
                graph->AP.child[x]->child[0]->checked = 1;
                graph->AP.child[x]->child[0]->pp[t] = graph->AP.child[x];
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->num_child > 0) {
            graph->AP.blockers[x]->num_child = 0;
            graph->AP.blockers[x]->child[0] = NULL;
        }
        graph->AP.blockers[x]->idx = 0;
        graph->AP.blockers[x]->num_parent = 0;
        memset(graph->AP.blockers[x]->parent, 0, sizeof(struct node *) * MAX_NODE);
        memset(graph->AP.blockers[x]->distance, 0, sizeof(double) * MAX_NODE);
        
        for (int y = 0; y < graph->AP.num_child; y++) { //Find parents and distance
            if (check_blockage_node(graph->AP.blockers[x], graph->AP.child[y], graph->coordinate, 0) == 0) {
                graph->AP.blockers[x]->parent[graph->AP.blockers[x]->num_parent] = graph->AP.child[y];
                graph->AP.blockers[x]->distance[graph->AP.blockers[x]->num_parent++] = distance(graph->AP.blockers[x], graph->AP.child[y]);
            }
        }
        
        struct node *temp_node = NULL;
        double temp;
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                for (int z = y + 1; z < graph->AP.blockers[x]->num_parent; z++) {
                    if (graph->AP.blockers[x]->parent[y]->height < graph->AP.blockers[x]->parent[z]->height) {
                        temp_node = graph->AP.blockers[x]->parent[y];
                        graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                        graph->AP.blockers[x]->parent[z] = temp_node;
                        temp = graph->AP.blockers[x]->distance[y];
                        graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                        graph->AP.blockers[x]->distance[z] = temp;
                    } else if (graph->AP.blockers[x]->parent[y]->height == graph->AP.blockers[x]->parent[z]->height) {
                        if (graph->AP.blockers[x]->distance[y] > graph->AP.blockers[x]->distance[z]) {
                            temp = graph->AP.blockers[x]->distance[y];
                            graph->AP.blockers[x]->distance[y] = graph->AP.blockers[x]->distance[z];
                            graph->AP.blockers[x]->distance[z] = temp;
                            temp_node = graph->AP.blockers[x]->parent[y];
                            graph->AP.blockers[x]->parent[y] = graph->AP.blockers[x]->parent[z];
                            graph->AP.blockers[x]->parent[z] = temp_node;
                        }
                    }
                }
            }
        }
    }
    
    return ret;
}

double greedy_matching(struct graph *graph, int t) {
    double z = 0;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        int y;
        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
            if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                graph->AP.blockers[x]->parent[y]->child[graph->AP.blockers[x]->parent[y]->num_child++] = graph->AP.blockers[x];
                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                z++;
                graph->AP.blockers[x]->checked = 1;
                break;
            }
        }
        
        if (y == graph->AP.blockers[x]->num_parent) {
            graph->AP.blockers[x]->reachability++;
        }
    }
    
    return (double) graph->AP.num_blockers - z;
}

double greedy_matching_depth2(struct graph *graph, double rem, int t) {
    double z = 0;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->checked == 0) {
            for (int y = 0; y < graph->AP.num_blockers; y++) {
                if (graph->AP.blockers[y]->checked == 1 &&
                    graph->AP.blockers[y]->num_child == 0) {
                    if (check_blockage_node(graph->AP.blockers[x], graph->AP.blockers[y], graph->coordinate, 0) == 0) {
                        graph->AP.blockers[y]->child[graph->AP.blockers[y]->num_child++] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[y];
                        z++;
                        break;
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
    }
    
    return rem - z;
}

double group_matching(struct graph *graph, int t) {
    double z = 0;
    int count = 0;
    while (count < graph->AP.num_blockers) {
        for (int i = 0; i < graph->AP.num_blockers; i++) {
            if (graph->AP.blockers[i]->checked == 0) {
                int j;
                for (j = 0; j < graph->AP.blockers[i]->num_blockers; j++) {
                    if (graph->AP.blockers[i]->blockers[j]->blocked == 0) {
                        if (check_blockage_node(graph->AP.blockers[i], graph->AP.blockers[i]->blockers[j], graph->coordinate, 0) == 0) {
                            if (graph->AP.blockers[i]->blockers[j]->num_child == 0) {
                                graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                graph->AP.blockers[i]->blockers[j]->num_child = 1;
                                graph->AP.blockers[i]->idx = j;
                                graph->AP.blockers[i]->checked = 1;
                                graph->AP.blockers[i]->traversed = 1;
                                graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                count++;
                                z++;
                                break;
                            } else {
                                if (graph->AP.blockers[i]->blockers[j]->child[0]->traversed == 1) {
                                    if (graph->AP.blockers[i]->distance_group[j] <
                                        graph->AP.blockers[i]->blockers[j]->child[0]->distance_group[graph->AP.blockers[i]->blockers[j]->child[0]->idx]) {
                                        graph->AP.blockers[i]->blockers[j]->child[0]->checked = 0;
                                        graph->AP.blockers[i]->blockers[j]->child[0]->traversed = 0;
                                        graph->AP.blockers[i]->blockers[j]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                        graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                        graph->AP.blockers[i]->checked = 1;
                                        graph->AP.blockers[i]->traversed = 1;
                                        graph->AP.blockers[i]->idx = j;
                                        break;
                                    } else if (graph->AP.blockers[i]->distance_group[j] ==
                                               graph->AP.blockers[i]->blockers[j]->child[0]->distance_group[graph->AP.blockers[i]->blockers[j]->child[0]->idx]) {
                                        if (graph->AP.blockers[i]->height > graph->AP.blockers[i]->blockers[j]->child[0]->height) {
                                            graph->AP.blockers[i]->blockers[j]->child[0]->checked = 0;
                                            graph->AP.blockers[i]->blockers[j]->child[0]->traversed = 0;
                                            graph->AP.blockers[i]->blockers[j]->child[0]->pp[t] = NULL;
                                            
                                            graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                            graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                            graph->AP.blockers[i]->checked = 1;
                                            graph->AP.blockers[i]->traversed = 1;
                                            graph->AP.blockers[i]->idx = j;
                                            break;
                                        }
                                    }
                                } else {
                                    graph->AP.blockers[i]->blockers[j]->child[0]->checked = 0;
                                    graph->AP.blockers[i]->blockers[j]->child[0]->traversed = 0;
                                    graph->AP.blockers[i]->blockers[j]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                    graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                    graph->AP.blockers[i]->checked = 1;
                                    graph->AP.blockers[i]->traversed = 1;
                                    graph->AP.blockers[i]->idx = j;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (j == graph->AP.blockers[i]->num_blockers) {
                    for (j = 0; j < graph->AP.blockers[i]->num_parent; j++) {
                        if (graph->AP.blockers[i]->parent[j]->num_child == 0) {
                            graph->AP.blockers[i]->parent[j]->child[0] = graph->AP.blockers[i];
                            graph->AP.blockers[i]->parent[j]->num_child = 1;
                            graph->AP.blockers[i]->checked = 1;
                            graph->AP.blockers[i]->idx = j;
                            graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->parent[j];
                            count++;
                            z++;
                            break;
                        } else if (graph->AP.blockers[i]->parent[j]->child[0]->traversed == 0) {
                            if (graph->AP.blockers[i]->distance[j] <
                                graph->AP.blockers[i]->parent[j]->child[0]->distance[graph->AP.blockers[i]->parent[j]->child[0]->idx]) {
                                graph->AP.blockers[i]->parent[j]->child[0]->checked = 0;
                                graph->AP.blockers[i]->parent[j]->child[0]->pp[t] = NULL;
                                graph->AP.blockers[i]->parent[j]->child[0] = graph->AP.blockers[i];
                                graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->parent[j];
                                graph->AP.blockers[i]->checked = 1;
                                graph->AP.blockers[i]->idx = j;
                                break;
                            } else if (graph->AP.blockers[i]->distance[j] ==
                                       graph->AP.blockers[i]->parent[j]->child[0]->distance[graph->AP.blockers[i]->parent[j]->child[0]->idx]) {
                                if (graph->AP.blockers[i]->height > graph->AP.blockers[i]->parent[j]->child[0]->height) {
                                    graph->AP.blockers[i]->parent[j]->child[0]->checked = 0;
                                    graph->AP.blockers[i]->parent[j]->child[0]->pp[t] = NULL;
                                    graph->AP.blockers[i]->parent[j]->child[0] = graph->AP.blockers[i];
                                    graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->parent[j];
                                    graph->AP.blockers[i]->checked = 1;
                                    graph->AP.blockers[i]->idx = j;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (j == graph->AP.blockers[i]->num_parent) {
                    count++;
                    graph->AP.blockers[i]->checked = 1;
                    graph->AP.blockers[i]->reachability++;
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double group_matching_fair(struct graph *graph, int t) {
    double z = 0;
    int count = 0;
    while (count < graph->AP.num_blockers) {
        for (int i = 0; i < graph->AP.num_blockers; i++) {
            if (graph->AP.blockers[i]->checked == 0) {
                int j;
                for (j = 0; j < graph->AP.blockers[i]->num_blockers; j++) {
                    if (graph->AP.blockers[i]->blockers[j]->blocked == 0) {
                        if (check_blockage_node(graph->AP.blockers[i], graph->AP.blockers[i]->blockers[j], graph->coordinate, 0) == 0) {
                            if (graph->AP.blockers[i]->blockers[j]->num_child == 0) {
                                graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                graph->AP.blockers[i]->blockers[j]->num_child = 1;
                                graph->AP.blockers[i]->idx = j;
                                graph->AP.blockers[i]->checked = 1;
                                graph->AP.blockers[i]->traversed = 1;
                                graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                count++;
                                z++;
                                break;
                            } else {
                                if (graph->AP.blockers[i]->blockers[j]->child[0]->traversed == 1) {
                                    if (graph->AP.blockers[i]->reachability > graph->AP.blockers[i]->blockers[j]->child[0]->reachability) {
                                        graph->AP.blockers[i]->blockers[j]->child[0]->checked = 0;
                                        graph->AP.blockers[i]->blockers[j]->child[0]->traversed = 0;
                                        graph->AP.blockers[i]->blockers[j]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                        graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                        graph->AP.blockers[i]->checked = 1;
                                        graph->AP.blockers[i]->traversed = 1;
                                        graph->AP.blockers[i]->idx = j;
                                        break;
                                    } else if (graph->AP.blockers[i]->reachability == graph->AP.blockers[i]->blockers[j]->child[0]->reachability) {
                                        if (graph->AP.blockers[i]->distance_group[j] <
                                            graph->AP.blockers[i]->blockers[j]->child[0]->distance_group[graph->AP.blockers[i]->blockers[j]->child[0]->idx]) {
                                            graph->AP.blockers[i]->blockers[j]->child[0]->checked = 0;
                                            graph->AP.blockers[i]->blockers[j]->child[0]->traversed = 0;
                                            graph->AP.blockers[i]->blockers[j]->child[0]->pp[t] = NULL;
                                            
                                            graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                            graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                            graph->AP.blockers[i]->checked = 1;
                                            graph->AP.blockers[i]->traversed = 1;
                                            graph->AP.blockers[i]->idx = j;
                                            break;
                                        } else if (graph->AP.blockers[i]->distance_group[j] ==
                                                   graph->AP.blockers[i]->blockers[j]->child[0]->distance_group[graph->AP.blockers[i]->blockers[j]->child[0]->idx]) {
                                            if (graph->AP.blockers[i]->height < graph->AP.blockers[i]->blockers[j]->child[0]->height) {
                                                graph->AP.blockers[i]->blockers[j]->child[0]->checked = 0;
                                                graph->AP.blockers[i]->blockers[j]->child[0]->traversed = 0;
                                                graph->AP.blockers[i]->blockers[j]->child[0]->pp[t] = NULL;
                                                
                                                graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                                graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                                graph->AP.blockers[i]->checked = 1;
                                                graph->AP.blockers[i]->traversed = 1;
                                                graph->AP.blockers[i]->idx = j;
                                                break;
                                            }
                                        }
                                    }
                                } else {
                                    graph->AP.blockers[i]->blockers[j]->child[0]->checked = 0;
                                    graph->AP.blockers[i]->blockers[j]->child[0]->traversed = 0;
                                    graph->AP.blockers[i]->blockers[j]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[i]->blockers[j]->child[0] = graph->AP.blockers[i];
                                    graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->blockers[j];
                                    graph->AP.blockers[i]->checked = 1;
                                    graph->AP.blockers[i]->traversed = 1;
                                    graph->AP.blockers[i]->idx = j;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (j == graph->AP.blockers[i]->num_blockers) {
                    for (j = 0; j < graph->AP.blockers[i]->num_parent; j++) {
                        if (graph->AP.blockers[i]->parent[j]->num_child == 0) {
                            graph->AP.blockers[i]->parent[j]->child[0] = graph->AP.blockers[i];
                            graph->AP.blockers[i]->parent[j]->num_child = 1;
                            graph->AP.blockers[i]->checked = 1;
                            graph->AP.blockers[i]->idx = j;
                            graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->parent[j];
                            count++;
                            z++;
                            break;
                        } else if (graph->AP.blockers[i]->parent[j]->child[0]->traversed == 0) {
                            if (graph->AP.blockers[i]->reachability > graph->AP.blockers[i]->parent[j]->child[0]->reachability) {
                                graph->AP.blockers[i]->parent[j]->child[0]->checked = 0;
                                graph->AP.blockers[i]->parent[j]->child[0]->pp[t] = NULL;
                                graph->AP.blockers[i]->parent[j]->child[0] = graph->AP.blockers[i];
                                graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->parent[j];
                                graph->AP.blockers[i]->checked = 1;
                                graph->AP.blockers[i]->idx = j;
                                break;
                            } else if (graph->AP.blockers[i]->reachability == graph->AP.blockers[i]->parent[j]->child[0]->reachability) {
                                if (graph->AP.blockers[i]->distance[j] <
                                    graph->AP.blockers[i]->parent[j]->child[0]->distance[graph->AP.blockers[i]->parent[j]->child[0]->idx]) {
                                    graph->AP.blockers[i]->parent[j]->child[0]->checked = 0;
                                    graph->AP.blockers[i]->parent[j]->child[0]->pp[t] = NULL;
                                    graph->AP.blockers[i]->parent[j]->child[0] = graph->AP.blockers[i];
                                    graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->parent[j];
                                    graph->AP.blockers[i]->checked = 1;
                                    graph->AP.blockers[i]->idx = j;
                                    break;
                                } else if (graph->AP.blockers[i]->distance[j] ==
                                           graph->AP.blockers[i]->parent[j]->child[0]->distance[graph->AP.blockers[i]->parent[j]->child[0]->idx]) {
                                    if (graph->AP.blockers[i]->height < graph->AP.blockers[i]->parent[j]->child[0]->height) {
                                        graph->AP.blockers[i]->parent[j]->child[0]->checked = 0;
                                        graph->AP.blockers[i]->parent[j]->child[0]->pp[t] = NULL;
                                        graph->AP.blockers[i]->parent[j]->child[0] = graph->AP.blockers[i];
                                        graph->AP.blockers[i]->pp[t] = graph->AP.blockers[i]->parent[j];
                                        graph->AP.blockers[i]->checked = 1;
                                        graph->AP.blockers[i]->idx = j;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (j == graph->AP.blockers[i]->num_parent) {
                    count++;
                    graph->AP.blockers[i]->checked = 1;
                    graph->AP.blockers[i]->reachability++;
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//checked = matched
double stable_matching(struct graph *graph, int t) {
    double z = 0;
    int count = 0;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->parent[y]->num_child = 1;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->idx = y;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                        count++;
                        z++;
                        break;
                    } else if (graph->AP.blockers[x]->distance[y] <
                               graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->idx = y;
                        break;
                    } else if (graph->AP.blockers[x]->distance[y] ==
                               graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                        if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->idx = y;
                            break;
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    count++;
                    graph->AP.blockers[x]->checked = 1;
                    graph->AP.blockers[x]->reachability++;
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double stable_matching_fair(struct graph *graph, int t) {
    double z = 0;
    int count = 0;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->parent[y]->num_child = 1;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->idx = y;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                        count++;
                        z++;
                        break;
                    } else if (graph->AP.blockers[x]->reachability > graph->AP.blockers[x]->parent[y]->child[0]->reachability) {
                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->idx = y;
                        break;
                    } else if (graph->AP.blockers[x]->reachability == graph->AP.blockers[x]->parent[y]->child[0]->reachability) {
                        if (graph->AP.blockers[x]->distance[y] <
                            graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->idx = y;
                            break;
                        } else if (graph->AP.blockers[x]->distance[y] ==
                                   graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                            if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->idx = y;
                                break;
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    count++;
                    graph->AP.blockers[x]->checked = 1;
                    graph->AP.blockers[x]->reachability++;
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//parent preference: height
//outdated
double stable_matching_height(struct graph *graph) {
    double z = 0;
    int count = 0;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->parent[y]->num_child++;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->idx = y;
                        count++;
                        z++;
                        break;
                    } else if (graph->AP.blockers[x]->height >
                               graph->AP.blockers[x]->parent[y]->child[0]->height) {
                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->idx = y;
                        break;
                    } else if (graph->AP.blockers[x]->height ==
                               graph->AP.blockers[x]->parent[y]->child[0]->height) {
                        if (graph->AP.blockers[x]->distance[y] <
                            graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->idx = y;
                            break;
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    count++;
                    graph->AP.blockers[x]->checked = 1;
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//outdated
double stable_matching_close(struct graph *graph) {
    double z = 0;
    int count = 0;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->height > graph->AP.blockers[x]->height) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child++;
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->checked = 1;
                            count++;
                            z++;
                            break;
                        } else if (graph->AP.blockers[x]->distance[y] <
                                   graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->idx = y;
                            break;
                        } else if (graph->AP.blockers[x]->distance[y] ==
                                   graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                            if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->idx = y;
                                break;
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    double height = graph->AP.blockers[x]->height;
                    int success = 0;
                    do {
                        height -= .1;
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            if (graph->AP.blockers[x]->parent[y]->height > height &&
                                graph->AP.blockers[x]->parent[y]->height <= height + 0.1) {
                                if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->parent[y]->num_child++;
                                    graph->AP.blockers[x]->idx = y;
                                    count++;
                                    z++;
                                    success = 1;
                                    break;
                                }
                            }
                        }
                        if (y < graph->AP.blockers[x]->num_parent) {
                            break;
                        }
                    } while (height >= MIN_HEIGHT);
                    if (success == 0) { //Worst case
                        /*
                        int y;
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            for (int temp = 0;
                                 temp < graph->AP.blockers[x]->parent[y]->child[0]->num_parent;
                                 temp++) {
                                if (graph->AP.blockers[x]->parent[y]->child[0]->parent[temp]->num_child == 0) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->parent[temp]->num_child++;
                                    graph->AP.blockers[x]->parent[y]->child[0]->parent[temp]->child[0] = graph->AP.blockers[x]->parent[y]->child[0];
                                    graph->AP.blockers[x]->parent[y]->child[0]->idx = temp;
                                    success = 1;
                                    break;
                                }
                            }
                            if (success == 1) {
                                break;
                            }
                        }
                        if (success == 1) {
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            z++;
                        }*/
                        count++;
                    }
                    graph->AP.blockers[x]->checked = 1;
                }
                
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
        graph->AP.blockers[x]->traversed = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//parent preference: height
//outdated
double stable_matching_close_height(struct graph *graph) {
    double z = 0;
    int count = 0;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->height > graph->AP.blockers[x]->height) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child++;
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->checked = 1;
                            count++;
                            z++;
                            break;
                        } else if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->idx = y;
                            break;
                        } else if (graph->AP.blockers[x]->height == graph->AP.blockers[x]->parent[y]->child[0]->height) {
                            if (graph->AP.blockers[x]->distance[y] <
                                graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->idx = y;
                                break;
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    double height = graph->AP.blockers[x]->height;
                    int success = 0;
                    do {
                        height -= .1;
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            if (graph->AP.blockers[x]->parent[y]->height > height &&
                                graph->AP.blockers[x]->parent[y]->height <= height + 0.1) {
                                if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->parent[y]->num_child++;
                                    graph->AP.blockers[x]->idx = y;
                                    count++;
                                    z++;
                                    success = 1;
                                    break;
                                }
                            }
                        }
                        if (y < graph->AP.blockers[x]->num_parent) {
                            break;
                        }
                    } while (height >= MIN_HEIGHT);
                    if (success == 0) { //Worst case
                        /*int y;
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            for (int temp = 0;
                                 temp < graph->AP.blockers[x]->parent[y]->child[0]->num_parent;
                                 temp++) {
                                if (graph->AP.blockers[x]->parent[y]->child[0]->parent[temp]->num_child == 0) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->parent[temp]->num_child++;
                                    graph->AP.blockers[x]->parent[y]->child[0]->parent[temp]->child[0] = graph->AP.blockers[x]->parent[y]->child[0];
                                    graph->AP.blockers[x]->parent[y]->child[0]->idx = temp;
                                    success = 1;
                                    break;
                                }
                            }
                            if (success == 1) {
                                break;
                            }
                        }
                        if (success == 1) {
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            z++;
                        }*/
                        count++;
                    }
                    graph->AP.blockers[x]->checked = 1;
                }
                
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
        graph->AP.blockers[x]->traversed = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_greedy(struct graph *graph, int t) {
    double z = 0.0;
    for (int x = 0; x < graph->AP.num_blockers; x++) { //If node blocked
        if (graph->AP.blockers[x]->checked == 0) { //If node hasn't been marked
            int y = 0;
            for (; y < graph->AP.blockers[x]->num_parent; y++) {
                if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                    graph->AP.blockers[x]->parent[y]->checked == 0) {//Get free parent
                    graph->AP.blockers[x]->parent[y]->child[graph->AP.blockers[x]->parent[y]->num_child++] = graph->AP.blockers[x];
                    z++;
                    //graph->AP.blockers[x]->marked = 1;
                    graph->AP.blockers[x]->blocked = 2;
                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                    break;
                }
            }
            
            if (y == graph->AP.blockers[x]->num_parent) { //If no free parent found, break previous link
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                        graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) {
                        //Modify child
                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                        
                        //Modify parent
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        
                        //Modify node
                        z++;
                        graph->AP.blockers[x]->blocked = 2;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                        break;
                    } else if (graph->AP.blockers[x]->parent[y]->checked == 1) {
                        //Modify Parent
                        graph->AP.blockers[x]->parent[y]->pp[t]->num_child = 0;
                        graph->AP.blockers[x]->parent[y]->pp[t]->child[0] = NULL;
                        
                        //Modify Child
                        graph->AP.blockers[x]->parent[y]->pp[t] = NULL;
                        graph->AP.blockers[x]->parent[y]->checked = 0;
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->parent[y]->num_child = 1;
                        
                        //Modify Node
                        z++;
                        graph->AP.blockers[x]->blocked = 2;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                        break;
                    }
                }
            }
            
            if (y == graph->AP.blockers[x]->num_parent) {
                graph->AP.blockers[x]->reachability++;
            }
        } else { // (2) Blocked and Marked node use the same parent to maintain stability
            z++;
            graph->AP.blockers[x]->blocked = 2;
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
        graph->people[x].marked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_greedy_stable(struct graph *graph, int t) {
    double z = 0.0;
    for (int x = 0; x < graph->AP.num_blockers; x++) { //If node blocked
        if (graph->AP.blockers[x]->checked == 0) { //If node hasn't been marked
            int y = 0;
            for (; y < graph->AP.blockers[x]->num_parent; y++) {
                if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                    graph->AP.blockers[x]->parent[y]->checked == 0) {//Get free parent
                    graph->AP.blockers[x]->parent[y]->child[graph->AP.blockers[x]->parent[y]->num_child++] = graph->AP.blockers[x];
                    z++;
                    graph->AP.blockers[x]->blocked = 2;
                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                    break;
                }
            }
            
            if (y == graph->AP.blockers[x]->num_parent) { //Find best link among breakable parent-child pair
                int found = 0;
                y = 0;
                for (int k = 0; k < graph->AP.blockers[x]->num_parent; k++) { //possible to add threshold on distance
                    if (graph->AP.blockers[x]->parent[k]->checked == 1) {
                        if (found == 0) {
                            found = 1;
                            y = k;
                        } else if (graph->AP.blockers[x]->parent[k]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                            y = k;
                        }
                    } else if (graph->AP.blockers[x]->parent[k]->num_child == 1 &&
                               graph->AP.blockers[x]->parent[k]->child[0]->blocked == 0) {
                        if (found == 0) {
                            found = 1;
                            y = k;
                        } else if (graph->AP.blockers[x]->parent[k]->child[0]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                            y = k;
                        }
                    }
                }
                
                if (found == 1) {
                    if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                        graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) {
                        //Modify child
                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                    } else if (graph->AP.blockers[x]->parent[y]->checked == 1) {
                        //Modify Parent
                        graph->AP.blockers[x]->parent[y]->pp[t]->num_child = 0;
                        graph->AP.blockers[x]->parent[y]->pp[t]->child[0] = NULL;
                        
                        //Modify Child
                        graph->AP.blockers[x]->parent[y]->pp[t] = NULL;
                        graph->AP.blockers[x]->parent[y]->checked = 0;
                        graph->AP.blockers[x]->parent[y]->num_child = 1;
                    }
                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                    z++;
                    graph->AP.blockers[x]->blocked = 2;
                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                } else {
                    graph->AP.blockers[x]->reachability++;
                }
            }
        } else { // (2) Blocked and Marked node use the same parent to maintain stability
            z++;
            graph->AP.blockers[x]->blocked = 2;
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//outdated
double update_depth2(struct graph *graph) {
    double z = 0.0;
    int count = 0;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->checked == 1) {
            z++;
            count++;
            graph->AP.blockers[x]->checked = 2;
            graph->AP.blockers[x]->blocked = 2;
        }
    }
    
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) { //If node blocked
            if (graph->AP.blockers[x]->checked == 0 && graph->AP.blockers[x]->marked == 0) { //If node hasn't been checked
                int y = 0, success = 0;
                for (; y < graph->AP.blockers[x]->num_parent; y++) { //For all parents in LOS of AP
                    if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                        graph->AP.blockers[x]->parent[y]->traversed == 0) {//Get free parent
                        graph->AP.blockers[x]->parent[y]->child[graph->AP.blockers[x]->parent[y]->num_child++] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->traversed = 1;
                        graph->AP.blockers[x]->checked = 1;
                        z++;
                        success = 1;
                        graph->AP.blockers[x]->blocked = 2;
                        break;
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {//Get Depth 1 parent
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[y]->traversed == 1) {
                            graph->AP.blockers[x]->parent[y]->child[graph->AP.blockers[x]->parent[y]->num_child++] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->traversed = 2;
                            z++;
                            graph->AP.blockers[x]->checked = 1;
                            success = 1;
                            graph->AP.blockers[x]->blocked = 2;
                            break;
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) { //No depth 1 LOS parent, find depth 1 blocked parent
                    for (y = 0; y < graph->AP.num_blockers; y++) {
                        if (graph->AP.blockers[y]->checked > 0 &&
                            graph->AP.blockers[y]->num_child == 0 &&
                            graph->AP.blockers[y]->traversed < 2) {
                            if (check_blockage_node(graph->AP.blockers[x], graph->AP.blockers[y], graph->coordinate, 0) == 0) {
                                graph->AP.blockers[y]->child[graph->AP.blockers[y]->num_child++] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->traversed = 2;
                                z++;
                                graph->AP.blockers[x]->checked = 1;
                                success = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                break;
                            }
                        }
                    }
                }
                
                if (success == 0) { //If no free parent found, kick tail node that has LOS
                    for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->num_child == 1 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0]->blocked == 0) {
                            
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0]->stability++;
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0]->traversed = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0]->checked = 0;
                            memset(graph->AP.blockers[x]->parent[y]->child[0]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0]->num_parent = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->traversed = 2;
                            z++;
                            graph->AP.blockers[x]->checked = 1;
                            success = 1;
                            graph->AP.blockers[x]->blocked = 2;
                            
                            for (int z = 0; z < graph->AP.num_blockers; z++) {
                                if (graph->AP.blockers[z]->checked == 0 && graph->AP.blockers[z]->marked == 1) {
                                    graph->AP.blockers[z]->marked = 0;
                                    count--;
                                }
                            }
                            break;
                        }
                    }
                }
                
                if (success == 0) { //Else, kick node that has los with blocked child not using previously established link
                    for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->num_child == 1 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0]->marked == 0 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0]->checked == 1) {
                            graph->AP.blockers[x]->child[0] = graph->AP.blockers[x]->parent[y]->child[0]->child[0];
                            graph->AP.blockers[x]->num_child = 1;
                            
                            graph->AP.blockers[x]->parent[y]->child[0]->stability++;
                            graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0] = NULL;
                            memset(graph->AP.blockers[x]->parent[y]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                            graph->AP.blockers[x]->parent[y]->child[0]->num_parent = 0;
                            
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->traversed = 1;
                            z++;
                            graph->AP.blockers[x]->checked = 1;
                            success = 1;
                            graph->AP.blockers[x]->blocked = 2;
                            
                            for (int z = 0; z < graph->AP.num_blockers; z++) {
                                if (graph->AP.blockers[z]->checked == 0 && graph->AP.blockers[z]->marked == 1) {
                                    graph->AP.blockers[z]->marked = 0;
                                    count--;
                                }
                            }
                            break;
                        }
                    }
                }
                
                if (success == 0) { //Else, kick node that has los with blocked child
                    for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->num_child == 1 &&
                            graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) {
                            graph->AP.blockers[x]->child[0] = graph->AP.blockers[x]->parent[y]->child[0]->child[0];
                            graph->AP.blockers[x]->child[0]->stability++;
                            graph->AP.blockers[x]->num_child = 1;
                            
                            graph->AP.blockers[x]->parent[y]->child[0]->stability++;
                            graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->child[0] = NULL;
                            memset(graph->AP.blockers[x]->parent[y]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                            graph->AP.blockers[x]->parent[y]->child[0]->num_parent = 0;
                            
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->traversed = 1;
                            z++;
                            graph->AP.blockers[x]->checked = 1;
                            success = 1;
                            graph->AP.blockers[x]->blocked = 2;
                            
                            for (int z = 0; z < graph->AP.num_blockers; z++) {
                                if (graph->AP.blockers[z]->checked == 0 && graph->AP.blockers[z]->marked == 1) {
                                    graph->AP.blockers[z]->marked = 0;
                                    count--;
                                }
                            }
                            break;
                        }
                    }
                }
                
                graph->AP.blockers[x]->marked = 1;
                count++;
            }
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].traversed = 0;
        graph->people[x].checked = 0;
        graph->people[x].marked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_group(struct graph *graph, double z, int t) {
    int count = ((int)(z + 0.5));
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->marked == 0 &&
                graph->AP.blockers[x]->checked == 0) {
                int success = 0;
                for (int i = 0; i < graph->AP.blockers[x]->num_blockers; i++) {
                    if (graph->AP.blockers[x]->blockers[i]->blocked == 0 &&
                        graph->AP.blockers[x]->blockers[i]->marked == 0) {
                        if (check_blockage_node(graph->AP.blockers[x], graph->AP.blockers[x]->blockers[i], graph->coordinate, 0) == 0) {
                            if (graph->AP.blockers[x]->blockers[i]->num_child == 0) {
                                graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->blockers[i]->num_child = 1;
                                graph->AP.blockers[x]->idx = i;
                                graph->AP.blockers[x]->traversed = 1;
                                graph->AP.blockers[x]->checked = 1;
                                count++;
                                z++;
                                success = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                break;
                            } else if (graph->AP.blockers[x]->blockers[i]->child[0]->marked == 0 &&
                                       graph->AP.blockers[x]->blockers[i]->child[0]->blocked > 0) {
                                if (graph->AP.blockers[x]->blockers[i]->child[0]->traversed == 1) {
                                    if (graph->AP.blockers[x]->distance_group[i] <
                                        graph->AP.blockers[x]->blockers[i]->child[0]->distance_group[graph->AP.blockers[x]->blockers[i]->child[0]->idx]) {
                                        graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->traversed = 1;
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->idx = i;
                                        success = 1;
                                        graph->AP.blockers[x]->blocked = 2;
                                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                        break;
                                    } else if (graph->AP.blockers[x]->distance_group[i] ==
                                               graph->AP.blockers[x]->blockers[i]->child[0]->distance_group[graph->AP.blockers[x]->blockers[i]->child[0]->idx]) {
                                        if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->blockers[i]->child[0]->height) {
                                            graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                            
                                            graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->traversed = 1;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->idx = i;
                                            success = 1;
                                            graph->AP.blockers[x]->blocked = 2;
                                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                            break;
                                        }
                                    }
                                } else {
                                    graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->idx = i;
                                    success = 1;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                    break;
                                }
                            }
                        }
                    }
                }
                //fprintf(stderr, "a\n");
                if (success == 0) {
                    for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->marked == 0) {
                            if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->parent[y]->num_child = 1;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 0;
                                graph->AP.blockers[x]->idx = y;
                                count++;
                                z++;
                                success = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                break;
                            } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 0 &&
                                       graph->AP.blockers[x]->parent[y]->child[0]->marked == 0) {
                                if (graph->AP.blockers[x]->distance[y] <
                                    graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 0;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                    success = 1;
                                    break;
                                } else if (graph->AP.blockers[x]->distance[y] ==
                                           graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->idx = y;
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->traversed = 0;
                                        graph->AP.blockers[x]->blocked = 2;
                                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                        success = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                //fprintf(stderr, "b\n");
                if (success == 0) {
                    for (int i = 0; i < graph->AP.blockers[x]->num_parent; i++) {
                        if (graph->AP.blockers[x]->parent[i]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[i]->marked == 1) {
                            graph->AP.blockers[x]->parent[i]->pp[t]->num_child = 0;
                            graph->AP.blockers[x]->parent[i]->pp[t]->child[0] = NULL;
                            graph->AP.blockers[x]->parent[i]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[i]->marked = 0;
                            graph->AP.blockers[x]->parent[i]->traversed = 0;
                            
                            graph->AP.blockers[x]->parent[i]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[i]->num_child = 1;
                            graph->AP.blockers[x]->idx = i;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 0;
                            count++;
                            z++;
                            success = 1;
                            graph->AP.blockers[x]->blocked = 2;
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[i];
                            break;
                        } else if (graph->AP.blockers[x]->parent[i]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[i]->child[0]->blocked == 0) {
                            graph->AP.blockers[x]->parent[i]->child[0]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[i]->child[0]->marked = 0;
                            graph->AP.blockers[x]->parent[i]->child[0]->traversed = 0;
                            
                            graph->AP.blockers[x]->parent[i]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->idx = i;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 0;
                            count++;
                            z++;
                            success = 1;
                            graph->AP.blockers[x]->blocked = 2;
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[i];
                            break;
                        }
                    }
                }
                //fprintf(stderr, "c\n");
                if (success == 0) {
                    graph->AP.blockers[x]->checked = 1;
                    count++;
                    graph->AP.blockers[x]->reachability++;
                }
            }
        }
        //fprintf(stderr, "count: %d, blockers: %d\n", count, graph->AP.num_blockers);
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
        graph->people[x].marked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_group_stable(struct graph *graph, double z, int t) {
    int count = ((int)(z + 0.5));
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->marked == 0 &&
                graph->AP.blockers[x]->checked == 0) {
                int success = 0;
                for (int i = 0; i < graph->AP.blockers[x]->num_blockers; i++) {
                    if (graph->AP.blockers[x]->blockers[i]->blocked == 0 &&
                        graph->AP.blockers[x]->blockers[i]->marked == 0) {
                        if (check_blockage_node(graph->AP.blockers[x], graph->AP.blockers[x]->blockers[i], graph->coordinate, 0) == 0) {
                            if (graph->AP.blockers[x]->blockers[i]->num_child == 0) {
                                graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->blockers[i]->num_child = 1;
                                graph->AP.blockers[x]->idx = i;
                                graph->AP.blockers[x]->traversed = 1;
                                graph->AP.blockers[x]->checked = 1;
                                count++;
                                z++;
                                success = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                break;
                            } else if (graph->AP.blockers[x]->blockers[i]->child[0]->marked == 0 &&
                                       graph->AP.blockers[x]->blockers[i]->child[0]->blocked > 0) {
                                if (graph->AP.blockers[x]->blockers[i]->child[0]->traversed == 1) {
                                    if (graph->AP.blockers[x]->distance_group[i] <
                                        graph->AP.blockers[x]->blockers[i]->child[0]->distance_group[graph->AP.blockers[x]->blockers[i]->child[0]->idx]) {
                                        graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->traversed = 1;
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->idx = i;
                                        success = 1;
                                        graph->AP.blockers[x]->blocked = 2;
                                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                        break;
                                    } else if (graph->AP.blockers[x]->distance_group[i] ==
                                               graph->AP.blockers[x]->blockers[i]->child[0]->distance_group[graph->AP.blockers[x]->blockers[i]->child[0]->idx]) {
                                        if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->blockers[i]->child[0]->height) {
                                            graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                            
                                            graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->traversed = 1;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->idx = i;
                                            success = 1;
                                            graph->AP.blockers[x]->blocked = 2;
                                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                            break;
                                        }
                                    }
                                } else {
                                    graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->idx = i;
                                    success = 1;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (success == 0) {
                    for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->marked == 0) {
                            if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->parent[y]->num_child = 1;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 0;
                                graph->AP.blockers[x]->idx = y;
                                count++;
                                z++;
                                success = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                break;
                            } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 0 &&
                                       graph->AP.blockers[x]->parent[y]->child[0]->marked == 0) {
                                if (graph->AP.blockers[x]->distance[y] <
                                    graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 0;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                    success = 1;
                                    break;
                                } else if (graph->AP.blockers[x]->distance[y] ==
                                           graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->idx = y;
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->traversed = 0;
                                        graph->AP.blockers[x]->blocked = 2;
                                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                        success = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (success == 0) {
                    int y = 0;
                    for (int k = 0; k < graph->AP.blockers[x]->num_parent; k++) { //possible to add threshold on distance
                        if (graph->AP.blockers[x]->parent[k]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[k]->marked == 1) {
                            if (success == 0) {
                                success = 1;
                                y = k;
                            } else if (graph->AP.blockers[x]->parent[k]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        } else if (graph->AP.blockers[x]->parent[k]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[k]->child[0]->blocked == 0) {
                            if (success == 0) {
                                success = 1;
                                y = k;
                            } else if (graph->AP.blockers[x]->parent[k]->child[0]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        }
                    }
                    
                    if (success == 1) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[y]->marked == 1) {
                            //Modify Parent
                            graph->AP.blockers[x]->parent[y]->pp[t]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t]->child[0] = NULL;
                            
                            //Modify Node
                            graph->AP.blockers[x]->parent[y]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->marked = 0;
                            graph->AP.blockers[x]->parent[y]->traversed = 0;
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) {
                            //Modify Child
                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->child[0]->marked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                        }
                        
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->idx = y;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->traversed = 0;
                        count++;
                        z++;
                        graph->AP.blockers[x]->blocked = 2;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                    } else {
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->reachability++;
                        count++;
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
        graph->people[x].marked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_group_fair(struct graph *graph, double z, int t) {
    int count = ((int)(z + 0.5));
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->marked == 0 &&
                graph->AP.blockers[x]->checked == 0) {
                int success = 0;
                for (int i = 0; i < graph->AP.blockers[x]->num_blockers; i++) {
                    if (graph->AP.blockers[x]->blockers[i]->blocked == 0 &&
                        graph->AP.blockers[x]->blockers[i]->marked == 0) {
                        if (check_blockage_node(graph->AP.blockers[x], graph->AP.blockers[x]->blockers[i], graph->coordinate, 0) == 0) {
                            if (graph->AP.blockers[x]->blockers[i]->num_child == 0) {
                                graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->blockers[i]->num_child = 1;
                                graph->AP.blockers[x]->idx = i;
                                graph->AP.blockers[x]->traversed = 1;
                                graph->AP.blockers[x]->checked = 1;
                                count++;
                                z++;
                                success = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                break;
                            } else if (graph->AP.blockers[x]->blockers[i]->child[0]->marked == 0 &&
                                       graph->AP.blockers[x]->blockers[i]->child[0]->blocked > 0) {
                                if (graph->AP.blockers[x]->blockers[i]->child[0]->traversed == 1) {
                                    
                                    if (graph->AP.blockers[x]->reachability > graph->AP.blockers[x]->blockers[i]->child[0]->reachability) {
                                        graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                        graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->traversed = 1;
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->idx = i;
                                        success = 1;
                                        graph->AP.blockers[x]->blocked = 2;
                                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                        break;
                                    } else if (graph->AP.blockers[x]->reachability == graph->AP.blockers[x]->blockers[i]->child[0]->reachability) {
                                        if (graph->AP.blockers[x]->distance_group[i] <
                                            graph->AP.blockers[x]->blockers[i]->child[0]->distance_group[graph->AP.blockers[x]->blockers[i]->child[0]->idx]) {
                                            graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                            graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                            
                                            graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->traversed = 1;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->idx = i;
                                            success = 1;
                                            graph->AP.blockers[x]->blocked = 2;
                                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                            break;
                                        } else if (graph->AP.blockers[x]->distance_group[i] ==
                                                   graph->AP.blockers[x]->blockers[i]->child[0]->distance_group[graph->AP.blockers[x]->blockers[i]->child[0]->idx]) {
                                            if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->blockers[i]->child[0]->height) {
                                                graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                                graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                                graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                                graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                                
                                                graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                                graph->AP.blockers[x]->traversed = 1;
                                                graph->AP.blockers[x]->checked = 1;
                                                graph->AP.blockers[x]->idx = i;
                                                success = 1;
                                                graph->AP.blockers[x]->blocked = 2;
                                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                                break;
                                            }
                                        }
                                    }
                                } else {
                                    graph->AP.blockers[x]->blockers[i]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->blockers[i]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->blockers[i]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->idx = i;
                                    success = 1;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->blockers[i];
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (success == 0) {
                    for (int y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->marked == 0) {
                            if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->parent[y]->num_child = 1;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 0;
                                graph->AP.blockers[x]->idx = y;
                                count++;
                                z++;
                                success = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                break;
                            } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 0 &&
                                       graph->AP.blockers[x]->parent[y]->child[0]->marked == 0) {
                                if (graph->AP.blockers[x]->reachability > graph->AP.blockers[x]->parent[y]->child[0]->reachability) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 0;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                    success = 1;
                                    break;
                                } else if (graph->AP.blockers[x]->reachability == graph->AP.blockers[x]->parent[y]->child[0]->reachability) {
                                    if (graph->AP.blockers[x]->distance[y] <
                                        graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->idx = y;
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->traversed = 0;
                                        graph->AP.blockers[x]->blocked = 2;
                                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                        success = 1;
                                        break;
                                    } else if (graph->AP.blockers[x]->distance[y] ==
                                               graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                        if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                            graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                            graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                            
                                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->idx = y;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->traversed = 0;
                                            graph->AP.blockers[x]->blocked = 2;
                                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                            success = 1;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (success == 0) {
                    int y = 0;
                    for (int k = 0; k < graph->AP.blockers[x]->num_parent; k++) { //possible to add threshold on distance
                        if (graph->AP.blockers[x]->parent[k]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[k]->marked == 1) {
                            if (success == 0) {
                                success = 1;
                                y = k;
                            } else if (graph->AP.blockers[x]->parent[k]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        } else if (graph->AP.blockers[x]->parent[k]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[k]->child[0]->blocked == 0) {
                            if (success == 0) {
                                success = 1;
                                y = k;
                            } else if (graph->AP.blockers[x]->parent[k]->child[0]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        }
                    }
                    
                    if (success == 1) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[y]->marked == 1) {
                            //Modify Parent
                            graph->AP.blockers[x]->parent[y]->pp[t]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t]->child[0] = NULL;
                            
                            //Modify Node
                            graph->AP.blockers[x]->parent[y]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->marked = 0;
                            graph->AP.blockers[x]->parent[y]->traversed = 0;
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) {
                            //Modify Child
                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->child[0]->marked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                        }
                        
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->idx = y;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->traversed = 0;
                        count++;
                        z++;
                        graph->AP.blockers[x]->blocked = 2;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                    } else {
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->reachability++;
                        count++;
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
        graph->people[x].marked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_stable(struct graph *graph, double z, int t) {
    double count = z;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {//Stable matching free links
                    if (graph->AP.blockers[x]->parent[y]->checked == 0) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            graph->AP.blockers[x]->idx = y;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                            break;
                        } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 1) {
                            if (graph->AP.blockers[x]->distance[y] <
                                graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->idx = y;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                break;
                            } else if (graph->AP.blockers[x]->distance[y] ==
                                       graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;

                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                    break;
                                }
                            }
                        }
                    }
                }
                //fprintf(stderr, "help\n");
                if (y == graph->AP.blockers[x]->num_parent) { //Find best link among breakable parent-child pair
                    for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                            //Modify Parent of Node using previous link
                            graph->AP.blockers[x]->parent[y]->pp[t]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t]->child[0] = NULL;
                            
                            //Modify Node using previous link
                            graph->AP.blockers[x]->parent[y]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t] = NULL;

                            //Update blocked node
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                            break;
                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                            //Modify child that is using previous link
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                            
                            //Update blocked node
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                            break;
                        }
                    }
                }
                //fprintf(stderr, "me\n");
                if (y == graph->AP.blockers[x]->num_parent) {
                    count++;
                    graph->AP.blockers[x]->checked = 1;
                }
            }
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
        graph->people[x].traversed = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_stable_stable(struct graph *graph, double z, int t) {
    double count = z;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {//Stable matching free links
                    if (graph->AP.blockers[x]->parent[y]->checked == 0) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            graph->AP.blockers[x]->idx = y;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                            break;
                        } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 1) {
                            if (graph->AP.blockers[x]->distance[y] <
                                graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->idx = y;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                break;
                            } else if (graph->AP.blockers[x]->distance[y] ==
                                       graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) { //Find best link among breakable parent-child pair
                    int found = 0;
                    y = 0;
                    for (int k = 0; k < graph->AP.blockers[x]->num_parent; k++) { //possible to add threshold on distance
                        if (graph->AP.blockers[x]->parent[k]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[k]->checked == 1) {
                            if (found == 0) {
                                found = 1;
                                y = k;
                            } else if (graph->AP.blockers[x]->parent[k]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        } else if (graph->AP.blockers[x]->parent[k]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[k]->child[0]->blocked == 0) {
                            if (found == 0) {
                                found = 1;
                                y = k;
                            } else if (graph->AP.blockers[x]->parent[k]->child[0]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        }
                    }
                    
                    if (found == 1) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                            //Modify Parent of Node using previous link
                            graph->AP.blockers[x]->parent[y]->pp[t]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t]->child[0] = NULL;
                            
                            //Modify Node using previous link
                            graph->AP.blockers[x]->parent[y]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                            //Modify child that is using previous link
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                        }
                        
                        //Update blocked node
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->idx = y;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->traversed = 1;
                        count++;
                        z++;
                        graph->AP.blockers[x]->blocked = 2;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                    } else {
                        count++;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->reachability++;
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
        graph->people[x].traversed = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_stable_fair(struct graph *graph, double z, int t) {
    double count = z;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {//Stable matching free links
                    if (graph->AP.blockers[x]->parent[y]->checked == 0) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            graph->AP.blockers[x]->idx = y;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                            break;
                        } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 1) {
                            if (graph->AP.blockers[x]->reachability > graph->AP.blockers[x]->parent[y]->child[0]->reachability) {
                                graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->idx = y;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                break;
                            } else if (graph->AP.blockers[x]->reachability == graph->AP.blockers[x]->parent[y]->child[0]->reachability) {
                                if (graph->AP.blockers[x]->distance[y] <
                                    graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->blocked = 2;
                                    graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                    break;
                                } else if (graph->AP.blockers[x]->distance[y] ==
                                           graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                        graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                                        
                                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->idx = y;
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->traversed = 1;
                                        graph->AP.blockers[x]->blocked = 2;
                                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) { //Find best link among breakable parent-child pair
                    int found = 0;
                    y = 0;
                    for (int k = 0; k < graph->AP.blockers[x]->num_parent; k++) { //possible to add threshold on distance
                        if (graph->AP.blockers[x]->parent[k]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[k]->checked == 1) {
                            if (found == 0) {
                                y = k;
                                found = 1;
                            } else if (graph->AP.blockers[x]->parent[k]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        } else if (graph->AP.blockers[x]->parent[k]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[k]->child[0]->blocked == 0) {
                            if (found == 0) {
                                y = k;
                                found = 1;
                            } else if (graph->AP.blockers[x]->parent[k]->child[0]->stability < graph->AP.blockers[x]->parent[y]->stability) {
                                y = k;
                            }
                        }
                    }
                    
                    if (found == 1) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                            //Modify Parent of Node using previous link
                            graph->AP.blockers[x]->parent[y]->pp[t]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t]->child[0] = NULL;
                            
                            //Modify Node using previous link
                            graph->AP.blockers[x]->parent[y]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->pp[t] = NULL;
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                            //Modify child that is using previous link
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->pp[t] = NULL;
                        }
                        
                        //Update blocked node
                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                        graph->AP.blockers[x]->idx = y;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->traversed = 1;
                        count++;
                        z++;
                        graph->AP.blockers[x]->blocked = 2;
                        graph->AP.blockers[x]->pp[t] = graph->AP.blockers[x]->parent[y];
                    } else {
                        count++;
                        graph->AP.blockers[x]->checked = 1;
                        graph->AP.blockers[x]->reachability++;
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
        graph->people[x].traversed = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//parent preference: height
//outdated
double update_stable_height(struct graph *graph, double z) {
    double count = z;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {//Stable matching free links
                    if (graph->AP.blockers[x]->parent[y]->checked == 0) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            break;
                        } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 1) {
                            if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->idx = y;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 1;
                                graph->AP.blockers[x]->blocked = 2;
                                break;
                            } else if (graph->AP.blockers[x]->height ==
                                       graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                if (graph->AP.blockers[x]->distance[y] <
                                    graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) { //Find best link among breakable parent-child pair
                    for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                            graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                            //Modify Parent of Node using previous link
                            graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->num_child = 0;
                            graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->child[0] = NULL;
                            
                            //Modify Node using previous link
                            graph->AP.blockers[x]->parent[y]->checked = 0;
                            memset(graph->AP.blockers[x]->parent[y]->parent, 0, sizeof(struct node *) * MAX_NODE);
                            memset(graph->AP.blockers[x]->parent[y]->distance, 0 , sizeof(double) * MAX_NODE);
                            graph->AP.blockers[x]->parent[y]->num_parent = 0;
                            graph->AP.blockers[x]->parent[y]->idx = 0;
                            graph->AP.blockers[x]->parent[y]->stability++;
                            
                            //Update blocked node
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            break;
                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                   graph->AP.blockers[x]->parent[y]->checked == 0 &&
                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                            //Modify child that is using previous link
                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                            memset(graph->AP.blockers[x]->parent[y]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                            memset(graph->AP.blockers[x]->parent[y]->child[0]->distance, 0, sizeof(double) * MAX_NODE);
                            graph->AP.blockers[x]->parent[y]->child[0]->num_parent = 0;
                            graph->AP.blockers[x]->parent[y]->child[0]->stability++;
                            //graph->AP.blockers[x]->parent[y]->num_child = 0;
                            //graph->AP.blockers[x]->parent[y]->child[0] = NULL;
                            
                            //Update blocked node
                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                            graph->AP.blockers[x]->idx = y;
                            graph->AP.blockers[x]->checked = 1;
                            graph->AP.blockers[x]->traversed = 1;
                            count++;
                            z++;
                            graph->AP.blockers[x]->blocked = 2;
                            break;
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    count++;
                    graph->AP.blockers[x]->checked = 1;
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
        graph->AP.blockers[x]->traversed = 0;
    }
    for (int x = 0; x < graph->AP.num_child; x++) {
        graph->AP.child[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//outdated
double update_stable_close(struct graph *graph, double z) {
    double count = z;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->height > graph->AP.blockers[x]->height) {
                        if (graph->AP.blockers[x]->parent[y]->checked == 0) {
                            if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->parent[y]->num_child = 1;
                                graph->AP.blockers[x]->idx = y;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 1;
                                count++;
                                z++;
                                graph->AP.blockers[x]->blocked = 2;
                                break;
                            } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 1) {
                                if (graph->AP.blockers[x]->distance[y] <
                                    graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                } else if (graph->AP.blockers[x]->distance[y] ==
                                           graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                    if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                        
                                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->traversed = 1;
                                        graph->AP.blockers[x]->idx = y;
                                        graph->AP.blockers[x]->blocked = 2;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    double height = graph->AP.blockers[x]->height;
                    do {
                        height -= .1;
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            if (graph->AP.blockers[x]->parent[y]->height > height &&
                                graph->AP.blockers[x]->parent[y]->height <= height + 0.1 &&
                                graph->AP.blockers[x]->parent[y]->checked == 0) {
                                if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->parent[y]->num_child = 1;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->idx = y;
                                    count++;
                                    z++;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                }
                            }
                        }
                        if (y < graph->AP.blockers[x]->num_parent) {
                            break;
                        }
                    } while (height >= MIN_HEIGHT);
                    
                    if (y == graph->AP.blockers[x]->num_parent) { //Worst case
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            if (graph->AP.blockers[x]->parent[y]->height > graph->AP.blockers[x]->height) {
                                if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                                    graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                                    //Modify Parent of Node using previous link
                                    graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->num_child = 0;
                                    graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->child[0] = NULL;
                                    
                                    //Modify Node using previous link
                                    graph->AP.blockers[x]->parent[y]->checked = 0;
                                    memset(graph->AP.blockers[x]->parent[y]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                    memset(graph->AP.blockers[x]->parent[y]->distance, 0 , sizeof(double) * MAX_NODE);
                                    graph->AP.blockers[x]->parent[y]->num_parent = 0;
                                    graph->AP.blockers[x]->parent[y]->idx = 0;
                                    graph->AP.blockers[x]->parent[y]->stability++;
                                    
                                    //Update blocked node
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->parent[y]->num_child = 1;
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    count++;
                                    z++;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                           graph->AP.blockers[x]->parent[y]->checked == 0 &&
                                           graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                                    //Modify child that is using previous link
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                    memset(graph->AP.blockers[x]->parent[y]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                    memset(graph->AP.blockers[x]->parent[y]->child[0]->distance, 0, sizeof(double) * MAX_NODE);
                                    graph->AP.blockers[x]->parent[y]->child[0]->num_parent = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->stability++;
                                    //graph->AP.blockers[x]->parent[y]->num_child = 0;
                                    //graph->AP.blockers[x]->parent[y]->child[0] = NULL;
                                    
                                    //Update blocked node
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    count++;
                                    z++;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                }
                            }
                        }
                        
                        if (y == graph->AP.blockers[x]->num_parent) {
                            double height = graph->AP.blockers[x]->height;
                            do {
                                height -= .1;
                                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                                    if (graph->AP.blockers[x]->parent[y]->height > height &&
                                        graph->AP.blockers[x]->parent[y]->height <= height + 0.1) {
                                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                                            graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                                            //Modify Parent of Node using previous link
                                            graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->num_child = 0;
                                            graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->child[0] = NULL;
                                            
                                            //Modify Node using previous link
                                            graph->AP.blockers[x]->parent[y]->checked = 0;
                                            memset(graph->AP.blockers[x]->parent[y]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                            memset(graph->AP.blockers[x]->parent[y]->distance, 0 , sizeof(double) * MAX_NODE);
                                            graph->AP.blockers[x]->parent[y]->num_parent = 0;
                                            graph->AP.blockers[x]->parent[y]->idx = 0;
                                            graph->AP.blockers[x]->parent[y]->stability++;
                                            
                                            //Update blocked node
                                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                                            graph->AP.blockers[x]->idx = y;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->traversed = 1;
                                            count++;
                                            z++;
                                            graph->AP.blockers[x]->blocked = 2;
                                            break;
                                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                                   graph->AP.blockers[x]->parent[y]->checked == 0 &&
                                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                                            //Modify child that is using previous link
                                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                            graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                            memset(graph->AP.blockers[x]->parent[y]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                            memset(graph->AP.blockers[x]->parent[y]->child[0]->distance, 0, sizeof(double) * MAX_NODE);
                                            graph->AP.blockers[x]->parent[y]->child[0]->num_parent = 0;
                                            graph->AP.blockers[x]->parent[y]->child[0]->stability++;
                                            //graph->AP.blockers[x]->parent[y]->num_child = 0;
                                            //graph->AP.blockers[x]->parent[y]->child[0] = NULL;
                                            
                                            //Update blocked node
                                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->idx = y;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->traversed = 1;
                                            count++;
                                            z++;
                                            graph->AP.blockers[x]->blocked = 2;
                                            break;
                                        }
                                    }
                                }
                                if (y < graph->AP.blockers[x]->num_parent) {
                                    break;
                                }
                            } while (height >= MIN_HEIGHT);
                        }
                        
                        if (y == graph->AP.blockers[x]->num_parent) {
                            count++;
                            graph->AP.blockers[x]->checked = 1;
                        }
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
        graph->AP.blockers[x]->traversed = 0;
    }
    for (int x = 0; x < graph->AP.num_child; x++) {
        graph->AP.child[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//parent preference: height
//outdated
double update_stable_close_height(struct graph *graph, double z) {
    double count = z;
    while (count < graph->AP.num_blockers) {
        for (int x = 0; x < graph->AP.num_blockers; x++) {
            if (graph->AP.blockers[x]->checked == 0) {
                int y;
                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                    if (graph->AP.blockers[x]->parent[y]->height > graph->AP.blockers[x]->height) {
                        if (graph->AP.blockers[x]->parent[y]->checked == 0) {
                            if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                graph->AP.blockers[x]->parent[y]->num_child = 1;
                                graph->AP.blockers[x]->idx = y;
                                graph->AP.blockers[x]->checked = 1;
                                graph->AP.blockers[x]->traversed = 1;
                                count++;
                                z++;
                                graph->AP.blockers[x]->blocked = 2;
                                break;
                            } else if (graph->AP.blockers[x]->parent[y]->child[0]->traversed == 1) {
                                if (graph->AP.blockers[x]->height > graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                    
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                } else if (graph->AP.blockers[x]->height ==
                                           graph->AP.blockers[x]->parent[y]->child[0]->height) {
                                    if (graph->AP.blockers[x]->distance[y] <
                                        graph->AP.blockers[x]->parent[y]->child[0]->distance[graph->AP.blockers[x]->parent[y]->child[0]->idx]) {
                                        graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->traversed = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                        graph->AP.blockers[x]->parent[y]->child[0]->blocked = 1;
                                        
                                        graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                        graph->AP.blockers[x]->checked = 1;
                                        graph->AP.blockers[x]->traversed = 1;
                                        graph->AP.blockers[x]->idx = y;
                                        graph->AP.blockers[x]->blocked = 2;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                
                if (y == graph->AP.blockers[x]->num_parent) {
                    double height = graph->AP.blockers[x]->height;
                    do {
                        height -= .1;
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            if (graph->AP.blockers[x]->parent[y]->height > height &&
                                graph->AP.blockers[x]->parent[y]->height <= height + 0.1 &&
                                graph->AP.blockers[x]->parent[y]->checked == 0) {
                                if (graph->AP.blockers[x]->parent[y]->num_child == 0) {
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->parent[y]->num_child = 1;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    graph->AP.blockers[x]->idx = y;
                                    count++;
                                    z++;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                }
                            }
                        }
                        if (y < graph->AP.blockers[x]->num_parent) {
                            break;
                        }
                    } while (height >= MIN_HEIGHT);
                    
                    if (y == graph->AP.blockers[x]->num_parent) { //Worst case
                        for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                            if (graph->AP.blockers[x]->parent[y]->height > graph->AP.blockers[x]->height) {
                                if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                                    graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                                    //Modify Parent of Node using previous link
                                    graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->num_child = 0;
                                    graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->child[0] = NULL;
                                    
                                    //Modify Node using previous link
                                    graph->AP.blockers[x]->parent[y]->checked = 0;
                                    memset(graph->AP.blockers[x]->parent[y]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                    memset(graph->AP.blockers[x]->parent[y]->distance, 0 , sizeof(double) * MAX_NODE);
                                    graph->AP.blockers[x]->parent[y]->num_parent = 0;
                                    graph->AP.blockers[x]->parent[y]->idx = 0;
                                    graph->AP.blockers[x]->parent[y]->stability++;
                                    
                                    //Update blocked node
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->parent[y]->num_child = 1;
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    count++;
                                    z++;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                           graph->AP.blockers[x]->parent[y]->checked == 0 &&
                                           graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                                    //Modify child that is using previous link
                                    graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                    memset(graph->AP.blockers[x]->parent[y]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                    memset(graph->AP.blockers[x]->parent[y]->child[0]->distance, 0, sizeof(double) * MAX_NODE);
                                    graph->AP.blockers[x]->parent[y]->child[0]->num_parent = 0;
                                    graph->AP.blockers[x]->parent[y]->child[0]->stability++;
                                    //graph->AP.blockers[x]->parent[y]->num_child = 0;
                                    //graph->AP.blockers[x]->parent[y]->child[0] = NULL;
                                    
                                    //Update blocked node
                                    graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                    graph->AP.blockers[x]->idx = y;
                                    graph->AP.blockers[x]->checked = 1;
                                    graph->AP.blockers[x]->traversed = 1;
                                    count++;
                                    z++;
                                    graph->AP.blockers[x]->blocked = 2;
                                    break;
                                }
                            }
                        }
                        
                        if (y == graph->AP.blockers[x]->num_parent) {
                            double height = graph->AP.blockers[x]->height;
                            do {
                                height -= .1;
                                for (y = 0; y < graph->AP.blockers[x]->num_parent; y++) {
                                    if (graph->AP.blockers[x]->parent[y]->height > height &&
                                        graph->AP.blockers[x]->parent[y]->height <= height + 0.1) {
                                        if (graph->AP.blockers[x]->parent[y]->num_child == 0 &&
                                            graph->AP.blockers[x]->parent[y]->checked == 1) { //Best link, force parent to use LOS
                                            //Modify Parent of Node using previous link
                                            graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->num_child = 0;
                                            graph->AP.blockers[x]->parent[y]->parent[graph->AP.blockers[x]->parent[y]->idx]->child[0] = NULL;
                                            
                                            //Modify Node using previous link
                                            graph->AP.blockers[x]->parent[y]->checked = 0;
                                            memset(graph->AP.blockers[x]->parent[y]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                            memset(graph->AP.blockers[x]->parent[y]->distance, 0 , sizeof(double) * MAX_NODE);
                                            graph->AP.blockers[x]->parent[y]->num_parent = 0;
                                            graph->AP.blockers[x]->parent[y]->idx = 0;
                                            graph->AP.blockers[x]->parent[y]->stability++;
                                            
                                            //Update blocked node
                                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->parent[y]->num_child = 1;
                                            graph->AP.blockers[x]->idx = y;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->traversed = 1;
                                            count++;
                                            z++;
                                            graph->AP.blockers[x]->blocked = 2;
                                            break;
                                        } else if (graph->AP.blockers[x]->parent[y]->num_child == 1 &&
                                                   graph->AP.blockers[x]->parent[y]->checked == 0 &&
                                                   graph->AP.blockers[x]->parent[y]->child[0]->blocked == 0) { //Best link, force parent to route to us
                                            //Modify child that is using previous link
                                            graph->AP.blockers[x]->parent[y]->child[0]->checked = 0;
                                            graph->AP.blockers[x]->parent[y]->child[0]->idx = 0;
                                            memset(graph->AP.blockers[x]->parent[y]->child[0]->parent, 0, sizeof(struct node *) * MAX_NODE);
                                            memset(graph->AP.blockers[x]->parent[y]->child[0]->distance, 0, sizeof(double) * MAX_NODE);
                                            graph->AP.blockers[x]->parent[y]->child[0]->num_parent = 0;
                                            graph->AP.blockers[x]->parent[y]->child[0]->stability++;
                                            //graph->AP.blockers[x]->parent[y]->num_child = 0;
                                            //graph->AP.blockers[x]->parent[y]->child[0] = NULL;
                                            
                                            //Update blocked node
                                            graph->AP.blockers[x]->parent[y]->child[0] = graph->AP.blockers[x];
                                            graph->AP.blockers[x]->idx = y;
                                            graph->AP.blockers[x]->checked = 1;
                                            graph->AP.blockers[x]->traversed = 1;
                                            count++;
                                            z++;
                                            graph->AP.blockers[x]->blocked = 2;
                                            break;
                                        }
                                    }
                                }
                                if (y < graph->AP.blockers[x]->num_parent) {
                                    break;
                                }
                            } while (height >= MIN_HEIGHT);
                        }
                        
                        if (y == graph->AP.blockers[x]->num_parent) {
                            count++;
                            graph->AP.blockers[x]->checked = 1;
                        }
                    }
                }
            }
        }
    }
    
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        graph->AP.blockers[x]->checked = 0;
        graph->AP.blockers[x]->traversed = 0;
    }
    for (int x = 0; x < graph->AP.num_child; x++) {
        graph->AP.child[x]->checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

void update_graph(int width, int length, struct graph *graph) {
    for (int x = 0; x < graph->population; x++) {
        int dir = rand() % 9;
        int increment = 1;
        switch (dir) {
            case 0:
                while (graph->people[x].y + increment < length) {
                    if (graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y + increment] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y + increment] = &graph->people[x];
                        graph->people[x].y += increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 1:
                while (graph->people[x].y + increment < length &&
                       graph->people[x].x + increment < width) {
                    if (graph->coordinate[(int)graph->people[x].x + increment][(int)graph->people[x].y + increment] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x + increment][(int)graph->people[x].y + increment] = &graph->people[x];
                        graph->people[x].y += increment;
                        graph->people[x].x += increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 2:
                while (graph->people[x].x + increment < width) {
                    if (graph->coordinate[(int)graph->people[x].x + increment][(int)graph->people[x].y] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x + increment][(int)graph->people[x].y] = &graph->people[x];
                        graph->people[x].x += increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 3:
                while (graph->people[x].y - increment >= 0 &&
                       graph->people[x].x + increment < width) {
                    if (graph->coordinate[(int)graph->people[x].x + increment][(int)graph->people[x].y - increment] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x + increment][(int)graph->people[x].y - increment] = &graph->people[x];
                        graph->people[x].y -= increment;
                        graph->people[x].x += increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 4:
                while (graph->people[x].y - increment >= 0) {
                    if (graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y - increment] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y - increment] = &graph->people[x];
                        graph->people[x].y -= increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 5:
                while (graph->people[x].y - increment >= 0 &&
                       graph->people[x].x - increment >= 0) {
                    if (graph->coordinate[(int)graph->people[x].x - increment][(int)graph->people[x].y - increment] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x - increment][(int)graph->people[x].y - increment] = &graph->people[x];
                        graph->people[x].y -= increment;
                        graph->people[x].x -= increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 6:
                while (graph->people[x].x - increment >= 0) {
                    if (graph->coordinate[(int)graph->people[x].x - increment][(int)graph->people[x].y] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x - increment][(int)graph->people[x].y] = &graph->people[x];
                        graph->people[x].x -= increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 7:
                while (graph->people[x].y + increment < length &&
                       graph->people[x].x - increment >= 0) {
                    if (graph->coordinate[(int)graph->people[x].x - increment][(int)graph->people[x].y + increment] == NULL) {
                        graph->coordinate[(int)graph->people[x].x][(int)graph->people[x].y] = NULL;
                        graph->coordinate[(int)graph->people[x].x - increment][(int)graph->people[x].y + increment] = &graph->people[x];
                        graph->people[x].y += increment;
                        graph->people[x].x -= increment;
                        break;
                    }
                    increment++;
                }
                break;
                
            case 8:
            default:
                break;
        }
    }
}

void update_graph_waypoint_group(int width, int length, struct graph *graph) {
    for (int x = 0; x < graph->population; x++) {
        if (graph->people[x].person == 0) {
            continue;
        }
        if ((graph->people[x].x != graph->people[x].x_dest) ||
            (graph->people[x].y != graph->people[x].y_dest)) {
            if (graph->people[x].x != graph->people[x].x_dest) {
                double increment = 1;
                if (graph->people[x].x_dest < graph->people[x].x) {
                    while (graph->people[x].x - increment >= 0) {
                        if (graph->coordinate[(int)(graph->people[x].x - increment + 0.5)][(int)(graph->people[x].y + 0.5)] == NULL) {
                            graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y + 0.5)] = NULL;
                            graph->coordinate[(int)(graph->people[x].x - increment + 0.5)][(int)(graph->people[x].y + 0.5)] = &graph->people[x];
                            graph->people[x].x -= increment;
                            break;
                        }
                        increment++;
                    }
                    if (graph->people[x].x - increment < 0) {
                        graph->people[x].timer++;
                    }
                } else {
                    while (graph->people[x].x + increment < width) {
                        if (graph->coordinate[(int)(graph->people[x].x + increment + 0.5)][(int)(graph->people[x].y + 0.5)] == NULL) {
                            graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y + 0.5)] = NULL;
                            graph->coordinate[(int)(graph->people[x].x + increment + 0.5)][(int)(graph->people[x].y + 0.5)] = &graph->people[x];
                            graph->people[x].x += increment;
                            break;
                        }
                        increment++;
                    }
                    if (graph->people[x].x + increment >= width) {
                        graph->people[x].timer++;
                    }
                }
            }
            if (graph->people[x].y != graph->people[x].y_dest) {
                double increment = 1;
                if (graph->people[x].y_dest < graph->people[x].y) {
                    while (graph->people[x].y - increment >= 0) {
                        if (graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y - increment + 0.5)] == NULL) {
                            graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y + 0.5)] = NULL;
                            graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y - increment + 0.5)] = &graph->people[x];
                            graph->people[x].y -= increment;
                            break;
                        }
                        increment++;
                    }
                    if (graph->people[x].y - increment < 0) {
                        graph->people[x].timer++;
                    }
                } else {
                    while (graph->people[x].y + increment < length) {
                        if (graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y + increment + 0.5)] == NULL) {
                            graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y + 0.5)] = NULL;
                            graph->coordinate[(int)(graph->people[x].x + 0.5)][(int)(graph->people[x].y + increment + 0.5)] = &graph->people[x];
                            graph->people[x].y += increment;
                            break;
                        }
                        increment++;
                    }
                    if (graph->people[x].y + increment >= length) {
                        graph->people[x].timer++;
                    }
                }
            }
            
            for (int y = 0; y < graph->people[x].num_blockers; y++) {
                int chx = 0, chy = 0;
                if (((int) (graph->people[x].x + graph->people[x].blockers[y]->x_dest + 0.5) >= 0) &&
                    ((int) (graph->people[x].x + graph->people[x].blockers[y]->x_dest + 0.5) < length)) {
                    if (graph->people[x].x + graph->people[x].blockers[y]->x_dest != graph->people[x].blockers[y]->x) {
                        chx = 1;
                    }
                }
                
                if (((int) (graph->people[x].y + graph->people[x].blockers[y]->y_dest + 0.5) >= 0) &&
                    ((int) (graph->people[x].y + graph->people[x].blockers[y]->y_dest + 0.5) < width)) {
                    if (graph->people[x].y + graph->people[x].blockers[y]->y_dest != graph->people[x].blockers[y]->y) {
                        chy = 1;
                    }
                }
                
                if (chx == 1) {
                    if (chy == 1) {
                        if (graph->coordinate[(int)(graph->people[x].x + graph->people[x].blockers[y]->x_dest + 0.5)][(int)(graph->people[x].y + graph->people[x].blockers[y]->y_dest + 0.5)] == NULL) {
                            graph->coordinate[(int)(graph->people[x].blockers[y]->x + 0.5)][(int)(graph->people[x].blockers[y]->y + 0.5)] = NULL;
                            graph->coordinate[(int)(graph->people[x].x + graph->people[x].blockers[y]->x_dest + 0.5)][(int)(graph->people[x].y + graph->people[x].blockers[y]->y_dest + 0.5)] = graph->people[x].blockers[y];
                            graph->people[x].blockers[y]->x = graph->people[x].x + graph->people[x].blockers[y]->x_dest;
                            graph->people[x].blockers[y]->y = graph->people[x].y + graph->people[x].blockers[y]->y_dest;
                        }
                    } else {
                        if (graph->coordinate[(int)(graph->people[x].x + graph->people[x].blockers[y]->x_dest + 0.5)][(int)(graph->people[x].blockers[y]->y + 0.5)] == NULL) {
                            graph->coordinate[(int)(graph->people[x].blockers[y]->x + 0.5)][(int)(graph->people[x].blockers[y]->y + 0.5)] = NULL;
                            graph->coordinate[(int)(graph->people[x].x + graph->people[x].blockers[y]->x_dest + 0.5)][(int)(graph->people[x].blockers[y]->y + 0.5)]= graph->people[x].blockers[y];
                            graph->people[x].blockers[y]->x = graph->people[x].x + graph->people[x].blockers[y]->x_dest;
                        }
                    }
                } else {
                    if (chy == 1) {
                        if (graph->coordinate[(int)(graph->people[x].blockers[y]->x + 0.5)][(int)(graph->people[x].y + graph->people[x].blockers[y]->y_dest + 0.5)] == NULL) {
                            graph->coordinate[(int)(graph->people[x].blockers[y]->x + 0.5)][(int)(graph->people[x].blockers[y]->y + 0.5)] = NULL;
                            graph->coordinate[(int)(graph->people[x].blockers[y]->x + 0.5)][(int)(graph->people[x].y + graph->people[x].blockers[y]->y_dest + 0.5)] = graph->people[x].blockers[y];
                            graph->people[x].blockers[y]->y = graph->people[x].y + graph->people[x].blockers[y]->y_dest;
                        }
                    }
                }
            }
            
            if (graph->people[x].timer >= 10) {
                graph->people[x].timer = rand() % 5;
                graph->people[x].x_dest = rand() % width;
                graph->people[x].y_dest = rand() % length;
            }
        } else if (graph->people[x].timer != 0) {
            if (graph->people[x].timer > 5) {
                graph->people[x].timer = rand() % 5;
            }
            graph->people[x].timer--;
        } else {
            graph->people[x].timer = rand() % 5;
            graph->people[x].x_dest = rand() % width;
            graph->people[x].y_dest = rand() % length;
        }
    }
}

void shift_index(struct graph *graph) {
    struct node *temp = graph->rr[0];
    for (int x = 0; x < graph->population - 1; x++) {
        graph->rr[x] = graph->rr[x+1];
    }
    graph->rr[graph->population - 1] = temp;
}

void sort_height_index(struct graph *graph) {
    struct node *temp = NULL;
    for (int x = 0; x < graph->population; x++) {
        for (int y = x + 1; y < graph->population; y++) {
            if (graph->rr[x]->height > graph->rr[y]->height) {
                temp = graph->rr[y];
                graph->rr[y] = graph->rr[x];
                graph->rr[x] = temp;
            }
        }
    }
}

void sort_stability(struct graph *graph) {
    struct node *temp = NULL;
    for (int x = 0; x < graph->population; x++) {
        for (int y = x + 1; y < graph->population; y++) {
            if (graph->rr[x]->stability < graph->rr[y]->stability) {
                temp = graph->rr[y];
                graph->rr[y] = graph->rr[x];
                graph->rr[x] = temp;
            }
        }
    }
}

void sort_reachability(struct graph *graph) {
    struct node *temp = NULL;
    for (int x = 0; x < graph->population; x++) {
        for (int y = x + 1; y < graph->population; y++) {
            if (graph->rr[x]->reachability < graph->rr[y]->reachability) {
                temp = graph->rr[y];
                graph->rr[y] = graph->rr[x];
                graph->rr[x] = temp;
            }
        }
    }
}

void reset(struct graph *graph) {
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].traversed = 0;
        graph->people[x].idx = 0;
        graph->people[x].checked = 0;
        graph->people[x].num_child = 0;
        graph->people[x].child[0] = NULL;
        memset(graph->people[x].parent, 0, sizeof(struct node *) * MAX_NODE);
        memset(graph->people[x].distance, 0, sizeof(double) * MAX_NODE);
        graph->people[x].num_parent = 0;
    }
}

int match(struct node *node, int t) {
    for (int x = 0; x < node->num_parent; x++) {
        if (node->parent[x]->traversed == 0) {
            node->parent[x]->traversed = 1;
            
            if (node->parent[x]->num_child == 0 || match(node->parent[x]->child[0], t)) {
                node->parent[x]->child[0] = node;
                node->pp[t] = node->parent[x];
                node->parent[x]->num_child = 1;
                return 1;
            }
        }
    }
    
    return 0;
}

int skip_match(struct node *node, int t) {
    for (int x = 0; x < node->num_parent; x++) {
        if ((node->parent[x]->num_child == 1 &&
             node->parent[x]->child[0]->checked == 1) ||
            (node->parent[x]->checked == 1)) {
            continue;
        }
        if (node->parent[x]->traversed == 0) {
            node->parent[x]->traversed = 1;
            
            if (node->parent[x]->num_child == 0 || skip_match(node->parent[x]->child[0], t)) {
                node->parent[x]->child[0] = node;
                node->parent[x]->num_child = 1;
                node->pp[t] = node->parent[x];
                node->blocked = 2;
                return 1;
            }
        }
    }
    
    return 0;
}

int all_match(struct node *node, int t) {
    if (node->blocked == 0 && node->checked == 1) {
        node->checked = 0;
        node->pp[t] = NULL;
        return 1;
    }
    
    for (int x = 0; x < node->num_parent; x++) {
        if (node->parent[x]->traversed == 0) {
            node->parent[x]->traversed = 1;
            
            if (node->parent[x]->num_child == 0 || all_match(node->parent[x]->child[0], t)) {
                node->parent[x]->pp[t] = NULL;
                node->parent[x]->child[0] = node;
                node->parent[x]->num_child = 1;
                node->pp[t] = node->parent[x];
                node->blocked = 2;
                return 1;
            }
        }
    }
    
    return 0;
}

int update_match(struct node *node, int t) {
    for (int x = 0; x < node->num_parent; x++) {
        if ((node->parent[x]->num_child == 1 &&
            node->parent[x]->child[0]->checked == 1) ||
            (node->parent[x]->checked == 1)) {
            continue;
        }
        if (node->parent[x]->traversed == 0) {
            node->parent[x]->traversed = 1;
            
            if (node->parent[x]->num_child == 0 || skip_match(node->parent[x]->child[0], t)) {
                node->parent[x]->child[0] = node;
                node->parent[x]->num_child = 1;
                node->pp[t] = node->parent[x];
                node->blocked = 2;
                return 1;
            }
        }
    }
    for (int x = 0; x < node->num_parent; x++) {
        if (node->parent[x]->traversed == 0) {
            node->parent[x]->traversed = 1;
            
            if (node->parent[x]->num_child == 0 || all_match(node->parent[x]->child[0], t)) {
                node->parent[x]->pp[t] = NULL;
                node->parent[x]->child[0] = node;
                node->parent[x]->num_child = 1;
                node->pp[t] = node->parent[x];
                node->blocked = 2;
                return 1;
            }
        }
    }

    return 0;
}

int update_match_stable(struct node *node, int t) {
    for (int x = 0; x < node->num_parent; x++) {
        if ((node->parent[x]->num_child == 1 &&
             node->parent[x]->child[0]->checked == 1) ||
            (node->parent[x]->checked == 1)) {
            continue;
        }
        if (node->parent[x]->traversed == 0) {
            node->parent[x]->traversed = 1;
            
            if (node->parent[x]->num_child == 0 || skip_match(node->parent[x]->child[0], t)) {
                node->parent[x]->child[0] = node;
                node->parent[x]->num_child = 1;
                node->pp[t] = node->parent[x];
                node->blocked = 2;
                return 1;
            }
        }
    }
    
    int found = 0, y = 0;
    for (int z = 0; z < node->num_parent; z++) { //possible to add threshold on distance
        if (node->parent[z]->checked == 1) {
            if (found == 0) {
                found = 1;
                y = z;
            } else if (node->parent[z]->stability < node->parent[y]->stability) {
                y = z;
            }
        } else if (node->parent[z]->num_child == 1 &&
                   node->parent[z]->child[0]->blocked == 0) {
            if (found == 0) {
                found = 1;
                y = z;
            } else if (node->parent[z]->child[0]->stability < node->parent[y]->stability) {
                y = z;
            }
        }
    }
    
    if (found == 1) {
        if (node->parent[y]->num_child == 1 &&
            node->parent[y]->child[0]->blocked == 0) {
            //Modify child
            node->parent[y]->child[0]->checked = 0;
            node->parent[y]->child[0]->pp[t] = NULL;
        } else if (node->parent[y]->checked == 1) {
            //Modify Parent
            node->parent[y]->pp[t]->num_child = 0;
            node->parent[y]->pp[t]->child[0] = NULL;
            
            //Modify Child
            node->parent[y]->pp[t] = NULL;
            node->parent[y]->checked = 0;
            node->parent[y]->num_child = 1;
        }
        node->parent[y]->child[0] = node;
        node->blocked = 2;
        node->pp[t] = node->parent[y];
        return 1;
    }
    
    return 0;
}

void reset_traversed(unsigned int pop, struct node *child[MAX_NODE]) {
    for (unsigned int i = 0; i < pop; i++) {
        child[i]->traversed = 0;
    }
}

double maximal_matching(struct graph *graph, int t) {
    double z = 0;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (match(graph->AP.blockers[x], t) == 1) {
            z++;
        } else {
            graph->AP.blockers[x]->reachability++;
        }
        reset_traversed(graph->AP.num_child, graph->AP.child);
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_perfect(struct graph *graph, int t) {
    double z = 0;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->checked == 0) {
            if (update_match(graph->AP.blockers[x], t) == 1) {
                z++;
                graph->AP.blockers[x]->blocked = 2;
            } else {
                graph->AP.blockers[x]->reachability++;
            }
            reset_traversed(graph->AP.num_child, graph->AP.child);
        } else {
            z++;
            graph->AP.blockers[x]->blocked = 2;
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

double update_perfect_stable(struct graph *graph, int t) {
    double z = 0;
    for (int x = 0; x < graph->AP.num_blockers; x++) {
        if (graph->AP.blockers[x]->checked == 0) {
            if (update_match_stable(graph->AP.blockers[x], t) == 1) {
                z++;
                graph->AP.blockers[x]->blocked = 2;
            } else {
                graph->AP.blockers[x]->reachability++;
            }
            reset_traversed(graph->AP.num_child, graph->AP.child);
        } else {
            z++;
            graph->AP.blockers[x]->blocked = 2;
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
    }
    
    return (double) graph->AP.num_blockers - z;
}

//used only for depth 2
double get_stability(struct graph *graph) {
    double ret = 0;
    for (int x = 0; x < graph->population; x++) {
        ret += graph->people[x].stability;
    }
    return ret;
}

double calc_stability(struct graph *graph, int t) {
    double ret = 0;
    for (int x = 0; x < graph->population; x++) {
        if (graph->people[x].pp[t] != graph->people[x].pp[t-1]) {
            graph->people[x].stability++;
        }
        ret += graph->people[x].stability;
    }
    return ret;
}

double get_jain(struct graph *graph) {
    double ret1 = 0;
    double ret2 = 0;
    for (int x = 0; x < graph->population; x++) {
        ret1 += graph->people[x].stability;
        ret2 += graph->people[x].stability * graph->people[x].stability;
    }
    ret1 = ret1 * ret1;
    ret2 = graph->population * ret2;
    if (ret2 != 0) {
        return ret1 / ret2;
    } else {
        return 1;
    }
}

double get_stabl(struct graph *graph) {
    double ret = 0;
    for (int x = 0; x < graph->population; x++) {
        if (graph->people[x].stability > ret) {
            ret = graph->people[x].stability;
        }
    }
    
    return ret;
}

double get_reach(struct graph *graph) {
    double ret = 0;
    for (int x = 0; x < graph->population; x++) {
        if (graph->people[x].reachability > ret) {
            ret = graph->people[x].reachability;
        }
    }
    
    return ret;
}

void update_capacity_delay(struct graph *graph, int t) {
    for (int x = 0; x < graph->population; x++) {
        if (graph->people[x].checked == 1) {
            continue;
        }
        if (graph->people[x].num_child > 0) {
            graph->people[x].child[0]->checked = 1;
            graph->people[x].delay = RENDER + NET + BEAM + 2 * IMAGE * 1000 / graph->people[x].capacity;
            graph->people[x].child[0]->capacity = calc_capacity(&graph->people[x], graph->people[x].child[0]);
            graph->people[x].child[0]->delay = graph->people[x].delay + BEAM;
            graph->people[x].capacity /= 2;
            if (graph->people[x].capacity > LIMIT) {
                graph->people[x].capacity = LIMIT;
                graph->people[x].child[0]->capacity = LIMIT;
            } else {
                graph->people[x].child[0]->capacity = graph->people[x].capacity;
            }
        } else if (graph->people[x].blocked == 0 && graph->people[x].pp[t] == NULL) {
            graph->people[x].delay = RENDER + NET + BEAM + IMAGE * 1000 / graph->people[x].capacity;
            if (graph->people[x].capacity > LIMIT) {
                graph->people[x].capacity = LIMIT;
            }
        }
    }
    
    for (int x = 0; x < graph->population; x++) {
        graph->people[x].checked = 0;
    }
}

double get_capacity(struct graph *graph) {
    double ret = 0;
    for (int x = 0; x < graph->population; x++) {
        ret += graph->people[x].capacity;
    }
    return ret/graph->population;
}

double get_delay(struct graph *graph) {
    double ret = 0;
    double count = 0;
    for (int x = 0; x < graph->population; x++) {
        if (graph->people[x].delay > 0) {
            ret += graph->people[x].delay;
            count++;
        }
    }
    return ret/count;
}

void init_stat(struct stat *stat) {
    for (int x = 0; x < 81; x++) {
        stat[x].stability = 0;
        stat[x].trials = 0;
    }
}

void save_stat(struct graph *graph, struct stat *stat) {
    for (int x = 0; x < graph->population; x++) {
        int idx = (int) (100 * (graph->people[x].height - MIN_HEIGHT) + 0.5);
        stat[idx].stability += graph->people[x].stability;
        stat[idx].trials++;
    }
}
















