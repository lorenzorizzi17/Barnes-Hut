#include "initial_distr.hpp"
#include "random"
#include "algorithm"


void makeCircularDistribution(std::vector<Particle>& particles, int N){
    // a central huge mass
    particles.push_back({ 0.1,  0.1,    0,   0,   1000000, 0, 0});
    // a goofy attempt to create a circular initial distribution of masses around a big stellar mass
    for(int i = 0; i < N-1; i++){
        double xcoord = 2 * std::cos(2 * M_PI/(double) N  * i) + 0.5 * (double)rand()/(double)RAND_MAX ;
        double ycoord = 2 * std::sin(2 * M_PI/(double) N  * i) + 0.5 * (double)rand()/(double)RAND_MAX ;
        double vx = -4 * std::sin(2 * M_PI/(double) N  * i);
        double vy = 4 * std::cos(2 * M_PI/(double) N  * i);
        particles.push_back({xcoord, ycoord, vx,vy, 2,0,0});
    } 
}

void makeRealisticGalaxy(std::vector<Particle>& particles, int N){
    static std::random_device rd;  
    std::mt19937 gen(rd()); 

    double lambda = 1/2.f;

    std::exponential_distribution<> exp_dist(lambda);

    double x_cent1 = 0.01;
    double y_cent1 = 0.01;
    particles.push_back({x_cent1, y_cent1, 0, 0, 1000000, 0, 0,0});
    double lim = 0.01;

    for (int i = 0; i < N-1; i++) {
        // Posizione radiale (distribuzione esponenziale, più densità al centro)
        double r = 1.5 + exp_dist(gen);// 1 + 4 * std::sqrt((double)rand()/RAND_MAX); 
        double theta = 2 * M_PI * (double)rand()/RAND_MAX;

        // Coordinate in x e y
        double xcoord = x_cent1 + r * std::cos(theta);
        double ycoord = y_cent1 + r * std::sin(theta);
        auto it = std::find_if(particles.begin(), particles.end(),[&](Particle const& p){return (std::abs(p.x-xcoord) < lim )&&(std::abs(p.y-ycoord) < lim); });
        while(it!= particles.end() || std::abs(xcoord) > 12 || std::abs(ycoord) > 12){
            r = 1.5 + exp_dist(gen);//1 + 8 * std::sqrt((double)rand() / RAND_MAX);
            theta = 2* M_PI * (double)rand() / RAND_MAX;
            xcoord = x_cent1 + r * std::cos(theta);
            ycoord = y_cent1 + r * std::sin(theta);
            it = std::find_if(particles.begin(), particles.end(),[&](Particle const& p){return (std::abs(p.x-xcoord) < lim )&&(std::abs(p.y-ycoord) < lim); });
        }

        double v= std::sqrt(G*1E6/r);
        double vx = -v * std::sin(theta);
        double vy = v * std::cos(theta);

        particles.push_back({xcoord, ycoord, vx, vy, 2, 0, 0, i+1});
    }
}

void makeDoubleGalaxies(std::vector<Particle>& particles, int N1, int N2){
    using namespace std;
    double x_cent1 = -4.8;
    double y_cent1 = 4.8;
    double x_cent2 = 4.8;
    double y_cent2 = -4.8;

    double lim = 0.05;

    // Corpo centrale per la prima galassia
    particles.push_back({x_cent1, y_cent1, 0.5, -1.1, 1000000, 0, 0, 0});

    // Distribuzione radiale e angolare delle stelle per la prima galassia
    for (int i = 0; i < N1 - 1; i++) {
        // Posizione radiale (distribuzione esponenziale, più densità al centro)
        double r = 1 + 3 * std::sqrt((double)rand() / RAND_MAX);  // raggio tra 0 e 2
        double theta = 2* M_PI * (double)rand() / RAND_MAX;  // Angolo casuale

        // Coordinate in x e y
        double xcoord = x_cent1 + r * std::cos(theta);
        double ycoord = y_cent1 + r * std::sin(theta);
        auto it = std::find_if(particles.begin(), particles.end(),[&](Particle const& p){return (std::abs(p.x-xcoord) < lim )&&(std::abs(p.y-ycoord) < lim); });
        while(it!= particles.end()){
            r = 1 + 3 * std::sqrt((double)rand() / RAND_MAX);
            theta = 2* M_PI * (double)rand() / RAND_MAX;
            xcoord = x_cent1 + r * std::cos(theta);
            ycoord = y_cent1 + r * std::sin(theta);
            it = std::find_if(particles.begin(), particles.end(),[&](Particle const& p){return (std::abs(p.x-xcoord) < lim )&&(std::abs(p.y-ycoord) < lim); });
        }

        double v= std::sqrt(G*1E6 / r);
        double vx = 0.5-v * std::sin(theta);
        double vy = -1.1+v * std::cos(theta);

        particles.push_back({xcoord, ycoord, vx, vy, 2, 0, 0, i+1});
    }

    particles.push_back({x_cent2, y_cent2, -0.5, 1.1, 1000000, 0, 0, N1 });

    // Distribuzione radiale e angolare delle stelle per la seconda galassia
    for (int i = 0; i < N2 - 1; i++) {
        // Posizione radiale (distribuzione esponenziale, più densità al centro)
        double r = 1 + 3 * std::sqrt((double)rand() / RAND_MAX);  // raggio tra 0 e 2
        double theta = 2* M_PI * (double)rand() / RAND_MAX;  // Angolo casuale

        // Coordinate in x e y
        double xcoord = x_cent2 + r * std::cos(theta);
        double ycoord = y_cent2 + r * std::sin(theta);
        auto it = std::find_if(particles.begin(), particles.end(),[&](Particle const& p){return (std::abs(p.x-xcoord) < lim )&&(std::abs(p.y-ycoord) < lim); });
        while(it!= particles.end()){
            r = 1 + 3* std::sqrt((double)rand() / RAND_MAX);
            theta = 2* M_PI * (double)rand() / RAND_MAX;
            xcoord = x_cent2 + r * std::cos(theta);
            ycoord = y_cent2 + r * std::sin(theta);
            it = std::find_if(particles.begin(), particles.end(),[&](Particle const& p){return (std::abs(p.x-xcoord) < lim )&&(std::abs(p.y-ycoord) < lim); });
        }

        double v= std::sqrt(G*1E6 / r);
        double vx = -0.5-v * std::sin(theta);
        double vy = 1.1+v * std::cos(theta);

        particles.push_back({xcoord, ycoord, vx, vy, 2, 0, 0, N1+i+1});
    }
}