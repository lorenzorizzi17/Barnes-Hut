#include"barnes-hutt.hpp"
#include"graphics.hpp"

#define NUM_THREADS 8


int main(){
    //main parameters
    double theta = 0.5;
    double dt = 0.001;
    const int N = 10000;    //10000
    omp_set_num_threads(NUM_THREADS);
    //graphical parameters
    int DIM = 12;

    //particles are collected in a std::vector 
    std::vector<Particle> particles;
    particles.reserve(N+1);

    //initial condition: create a circular-ish distribution
    makeRealisticGalaxy(particles, N);

    int time = 0;
    double start = omp_get_wtime();
    //graphical window and main loop
    while (time < 1000) {

        // -- dynamics --
        barnesHuttAlgorithm(particles, theta, dt, time);
        //nestedLoopAlgorithm(particles,dt,time );


        // -- remove particles that crossed the borders --
        for(auto it = particles.begin(); it != particles.end();){
            if(it->x <= -DIM || it->x >= DIM || it->y >= DIM || it->y <= -DIM){
                it = particles.erase(it);
            } else {
                it++;
            }
        }
        time++;
    }
        
    std::cout << omp_get_wtime()-start << "\n";
}

