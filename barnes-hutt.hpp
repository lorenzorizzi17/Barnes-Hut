#include<iostream>
#include<list>
#include<cassert>
#include<stack>
#include<list>
#include<math.h>
#include<omp.h>
#include<vector>


struct Particle{
    double x;
    double y;
    double v_x;
    double v_y;
    double m;
    double a_x = 0;
    double a_y = 0;
    int id;
};

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

void makeCircularDistribution(std::vector<Particle>&, int);

void makeDoubleGalaxies(std::vector<Particle>& particles, int, int);

void makeRealisticGalaxy(std::vector<Particle>&, int);

void init_node(node* & n , Particle* P, double xmin, double xmax, double ymin, double ymax);

void print_all_M(node* ROOT);

void adjust(int a, double & x_min, double & x_max, double & y_min, double & y_max);

node* buildQuadTree(std::vector<Particle>&);

node* buildQuadTreeParallel(std::list<Particle>&);

void freeQuadTree(node* ROOT);

double width(node const*const n);

double distance (Particle const part_1, double x_cm, double y_cm);

void force (double* forc, Particle const part_1, double mass, double x_cm, double y_cm);

void calculate_force(Particle& part, node*& n, double theta);

void calculate_force(Particle&, Particle&);

void leap_frog_step(Particle& part, double dt);

void leap_frog(std::vector<Particle>& particles, node* ROOT, double dt, int T);

void barnesHuttAlgorithm(std::vector<Particle>&, double, double,int);

void nestedLoopAlgorithm(std::vector<Particle>&,double,int);