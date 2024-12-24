#ifndef NBODY
#define NBODY

#include"gravity.hpp"

#define NUM_THREADS 8

enum simulationMode {BARNESHUTT, DOUBLELOOP};

enum initialDistribution {SINGLEGALAXY, DOUBLEGALAXY};

class Simulation{
    private:
    std::vector<Particle> m_particles;
    int m_time;

    public:
    void init(initialDistribution, int);
    void run(simulationMode, double);
};

#endif