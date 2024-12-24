//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//// The actual implementation of the Barnes-Hut algorithm, with the creation ////
////        of the quad-tree and later calculation of the forces              ////
////          Those tasks can be parallelized using OpenMP                    ////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

#ifndef BH
#define BH

#include "gravity.hpp"

struct node {
    double M;
    double x_cm;
    double y_cm;
    node* children [4] = {NULL,NULL,NULL,NULL};
    Particle* p;
    bool hasParticles = false;
    double x_min, x_max;
    double y_min, y_max;
    std::vector<Particle*> part; 
};

void barnesHuttAlgorithm(std::vector<Particle>&, double, double, int);

void initNode(node* & n , Particle* P, double xmin, double xmax, double ymin, double ymax);

void print_all_M(node* ROOT);

void adjust(int a, double & x_min, double & x_max, double & y_min, double & y_max);

node* buildQuadTree(std::vector<Particle>&);

node* buildQuadTreeParallel(std::vector<Particle>&);

void freeQuadTree(node* ROOT);

double width(node const*const n);

int getQuadrant(Particle, node*);

void calculate_force(Particle& part, node*& n, double theta);


#endif