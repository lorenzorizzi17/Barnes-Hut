//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//// Setting the basic parameter for the n-body problem and the naive way to  ////
////                    solve it using a dumb double loop                     ////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


#ifndef GRAVITY
#define GRAVITY

#include<iostream>
#include<cassert>
#include<stack>
#include<math.h>
#include<omp.h>
#include<vector>

#define G 3.4E-5

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


void calculate_force(Particle&, Particle&);

void eulerStep(Particle& part, double dt);

double distance(Particle const part_1, double x_cm, double y_cm);

void force(double* forc, Particle const part_1, double mass, double x_cm, double y_cm);

void nestedLoopAlgorithm(std::vector<Particle>&, double, int);



#endif

