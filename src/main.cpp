#include"nbody-sim.hpp"

int main(){
    //main parameters
    double theta = 0.5;
    double dt = 0.005;
    const int N = 20000;

    Simulation sim;
    sim.init(SINGLEGALAXY, N);
    sim.run(BARNESHUTT, dt);
}

