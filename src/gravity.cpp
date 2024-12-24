#include"gravity.hpp"

void calculate_force(Particle& p1, Particle& p2){
    double forc[2];
    force(forc, p1, p2.m, p2.x, p2.y);
    p1.a_x += forc[0]/p1.m; 
    p1.a_y += forc[1]/p1.m; 
    p2.a_x += - forc[0]/p2.m; 
    p2.a_y += - forc[1]/p2.m; 
}

double distance (Particle const part_1, double x_cm, double y_cm){
    return std::sqrt( std::pow(part_1.x - x_cm, 2) + std::pow(part_1.y - y_cm, 2) );
}

void force(double* forc, Particle const part_1, double mass, double x_cm, double y_cm){
    double r = distance(part_1, x_cm, y_cm); 
    assert(r != 0);
    *forc = G * (part_1.m * mass) / (r*r*r) * (x_cm - part_1.x);
    *(forc+1) = G * (part_1.m * mass)/(r*r*r) * (y_cm - part_1.y);
}

void eulerStep(Particle& part, double dt){
        part.x = part.x + part.v_x * dt; 
        part.y = part.y + part.v_y * dt; 
        part.v_x = part.v_x + part.a_x * dt; 
        part.v_y = part.v_y + part.a_y * dt; 
}

void nestedLoopAlgorithm(std::vector<Particle>& particles, double dt, int t){
    //first iteration of the leap-frog iteration
    if(t == 0){
         for (auto it = particles.begin(); it != particles.end(); it++) {
            for (auto it2 = it; it2 != particles.end(); it2++) {
                if (it != it2) {
                    calculate_force(*it, *it2);
                }
            }
        }
        for (auto it = particles.begin(); it != particles.end(); it++) {
            it->v_x += it->a_x * dt / 2.;
            it->v_y += it->a_y * dt / 2.;
            it->x += it->v_x * dt / 2.;
            it->y += it->v_y * dt / 2.;
        }
    }

    //calculate the forces
    #pragma omp parallel for
    for(int i = 0; i < particles.size(); i++){
        for(int j = i+1; j < particles.size(); j++){
            calculate_force(particles[i], particles[j]); 
        }
    }

    //update the bodies positions
    for(auto it = particles.begin(); it != particles.end(); it++){
        eulerStep(*it, dt);
    }
}