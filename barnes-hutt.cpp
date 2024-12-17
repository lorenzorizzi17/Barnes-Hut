#include"barnes-hutt.hpp"
#include<cstdlib> 
#include<algorithm>
#include <random>

#define G 3.4E-5


void makeCircularDistribution(std::vector<Particle>& particles, int N){
    // a central huge mass
    particles.push_back({ 0.1,  0.1,    0,   0,   1000000, 0, 0});
    // a goofy attempt to create a circular initial distribution of masses around a big stellar mass
    for(int i = 0; i < N-1; i++){
        double xcoord = 2*std::cos(2 * M_PI/(double) N  * i) + 0.5 * (double)rand()/(double)RAND_MAX ;
        double ycoord = 2*std::sin(2 * M_PI/(double) N  * i) + 0.5 * (double)rand()/(double)RAND_MAX ;
        double vx = -4*std::sin(2 * M_PI/(double) N  * i);
        double vy = 4*std::cos(2 * M_PI/(double) N  * i);
        particles.push_back({xcoord, ycoord, vx,vy, 2,0,0});
    } 
}

void makeRealisticGalaxy(std::vector<Particle>& particles, int N){
    static std::random_device rd;  
    std::mt19937 gen(rd());  // Generatore di numeri casuali

    // Imposta il parametro della distribuzione esponenziale (lambda)
    double lambda = 1/2.f;

    // Crea un oggetto per la distribuzione esponenziale con parametro lambda
    std::exponential_distribution<> exp_dist(lambda);

    double x_cent1 = 0.01;
    double y_cent1 = 0.01;
    particles.push_back({x_cent1, y_cent1, 0, 0, 1000000, 0, 0,0});
    double lim = 0.01;

    for (int i = 0; i < N-1; i++) {
        // Posizione radiale (distribuzione esponenziale, più densità al centro)
        double r = 0.5+ exp_dist(gen);// 1 + 4 * std::sqrt((double)rand()/RAND_MAX); 
        double theta = 2 * M_PI * (double)rand()/RAND_MAX;

        // Coordinate in x e y
        double xcoord = x_cent1 + r * std::cos(theta);
        double ycoord = y_cent1 + r * std::sin(theta);
        auto it = std::find_if(particles.begin(), particles.end(),[&](Particle const& p){return (std::abs(p.x-xcoord) < lim )&&(std::abs(p.y-ycoord) < lim); });
        while(it!= particles.end() || std::abs(xcoord) > 12 || std::abs(ycoord) > 12){
            r = 0.5 + exp_dist(gen);//1 + 8 * std::sqrt((double)rand() / RAND_MAX);
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

int getQuadrant(Particle p, node* curr_cell){
    double x = p.x;
    double y = p.y;
    double x_min = curr_cell->x_min;
    double x_max = curr_cell->x_max;
    double y_max = curr_cell->y_max;
    double y_min = curr_cell->y_min;
    double center_x = (x_min + x_max) / 2;
    double center_y = (y_min + y_max) / 2;
    assert(!std::isnan(x));
    if (x >= center_x && y >= center_y) {
        return 0; 
    } else if (x >= center_x && y <= center_y) {
        return 1; 
    } else if (x <= center_x && y <= center_y) {
        return 2; 
    } else if (x <= center_x && y >= center_y) {
        return 3;  
    }
    throw std::runtime_error{"getQuadrant threw an exception"};
};

//initialize a node with a given particle
void init_node(node* & n , Particle* P, double xmin, double xmax, double ymin, double ymax){
    n = new node;
    n->M = P->m;
    n->p = P;    
    n->x_cm = P->x;
    n->y_cm = P->y;
    for (int i = 0; i < 4; i++) {
        (n->children)[i] = nullptr;
    }  
    n->hasParticles = true;
    n->x_max = xmax;
    n->y_max = ymax;
    n->x_min = xmin;
    n->y_min = ymin;
}


void print_all_M(node* ROOT) {
    if (ROOT == nullptr) {
        return;  
    }
    
    std::stack<node*> toVisit;
    toVisit.push(ROOT);
    
    while (!toVisit.empty()) {
        node* current = toVisit.top();
        toVisit.pop();
        
        std::cout << "ROOT, content:" << (current->part).size() << ", " << current->hasParticles <<" , M: " << current->M << " CM (" << current->x_cm << ", " << current->y_cm << "), x: " << current->x_min << " - " <<  current->x_max << ", y: " <<  current->y_min << " - "  << current->y_max << std::endl;
        for (int i = 0; i < 4; ++i) {
            if (current->children[i] != nullptr) {
                toVisit.push(current->children[i]);
            }
        }
    }
}


void adjust(int a, double & x_min, double & x_max, double & y_min, double & y_max){
    switch(a){
        case 0:
            x_min = x_min + 0.5*(x_max-x_min);
            y_min = y_min + 0.5*(y_max-y_min);
            break;
        case 1:
            x_min = x_min + 0.5*(x_max-x_min);
            y_max = y_max - 0.5*(y_max-y_min);
            break;
        case 2:
            x_max = x_max - 0.5*(x_max-x_min);
            y_max = y_max - 0.5*(y_max-y_min);
            break;
        case 3:
            x_max = x_max - 0.5*(x_max-x_min);
            y_min = y_min + 0.5*(y_max-y_min);
            break;
        default: 
            break;
        }
}

//this function splits the bodies into the four children and launch 4 tasks for every children 
void split_part(node*& curr_cell){
  
    int N_part = curr_cell->part.size();

    if(N_part > 1){
        curr_cell->hasParticles = false;
    }

    for(int j = 0; j <N_part; j++){
        //starting from the root node
        Particle* p = (curr_cell->part[j]);
        double m = p->m;
        double x = p->x;
        double y = p->y;

        int a = getQuadrant(*p, curr_cell); 

        if(!curr_cell->children[a]){
            switch(a){
                case 0:
                    init_node(curr_cell->children[a], p, curr_cell->x_min + 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min + 0.5*(curr_cell->y_max-curr_cell->y_min), curr_cell->y_max ); 
                    break;
                case 1:
                    init_node(curr_cell->children[a], p, curr_cell->x_min +0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                    break;
                case 2:
                    init_node(curr_cell->children[a], p, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                    break;
                case 3:
                    init_node(curr_cell->children[a], p, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min + 0.5*(curr_cell->y_max -curr_cell->y_min ), curr_cell->y_max ); 
                    break;
                default: 
                    break;
                }
            } else {
                curr_cell->children[a]->x_cm = (curr_cell->children[a]->x_cm*curr_cell->children[a]->M + m*x) /(m+curr_cell->children[a]->M);
                curr_cell->children[a]->y_cm = (curr_cell->children[a]->y_cm*curr_cell->children[a]->M + m*y) /(m+curr_cell->children[a]->M);
                curr_cell->children[a]->M = curr_cell->children[a]->M + m; 
            }

        curr_cell->children[a]->part.push_back(p);
    }

    
    for (int i = 0; i < 4; i++) {
        if (curr_cell->children[i] && curr_cell->children[i]->part.size() > 1) {
            if(curr_cell->children[i]->part.size() < 300){
                split_part(curr_cell->children[i]);
            } else {
                #pragma omp task
                    split_part(curr_cell->children[i]);
            }
        }
    }

    (curr_cell->part).clear();
}

node* buildQuadTreeParallel(std::vector<Particle>& particles){

    node* ROOT = NULL;
    init_node(ROOT, &*particles.begin(), -12, 12, -12, 12);
    ROOT->part.push_back(&*particles.begin());

    for (auto it = particles.begin()+1; it != particles.end(); it++) {
        ROOT->part.push_back(&*it);
        ROOT->x_cm = (ROOT->x_cm*ROOT->M + (it)->m * (it)->x) /((it)->m+ROOT->M);
        ROOT->y_cm = (ROOT->y_cm*ROOT->M + (it)->m*(it)->y) /((it)->m+ROOT->M);
        ROOT->M = ROOT->M + (it)->m;
    }

    #pragma omp parallel
    {
        #pragma omp single
            split_part(ROOT);  
    }
    
    return ROOT;
}

node* buildQuadTree(std::vector<Particle>& particles){
    node* ROOT = NULL;
    //initialize the first cell with the fisrt particle
    init_node(ROOT, &particles.front(), -12, 12, -12, 12);

    //now run the algorithm for every particles
    auto it = particles.begin();
    it++;
    for(; it != particles.end(); it++){
        //starting from the root node
        node* curr_cell = ROOT;
        double m = it->m;
        double x = it->x;
        double y = it->y;

        double x_min = -12;
        double x_max = 12;
        double y_min = -12;
        double y_max = 12;

        //enter the loop until a break statement is found
        while(true){
            int b;
            if(curr_cell->hasParticles){    // enter this branch if the current node has already a particle (hence it's an external node)

                //the external node will become an internal one. Compute the total mass and CM
                curr_cell->x_cm = (curr_cell->x_cm*curr_cell->M + m*x) /(m+curr_cell->M);
                curr_cell->y_cm = (curr_cell->y_cm*curr_cell->M + m*y) /(m+curr_cell->M);
                curr_cell->M = curr_cell->M + m;

                //we need to first move the old particle
                Particle* old_particle = (curr_cell->p);
                // and see in which quadrant it belongs. a is just an int
                int a = getQuadrant(*old_particle, curr_cell);
                
                //the current cell becomes an internal node, so set hasParticle as false
                curr_cell->p = nullptr;
                curr_cell->hasParticles=false;

                //now depending on a, create a new node and put the old particle there
                switch(a){
                    case 0:
                        init_node(curr_cell->children[a], old_particle, curr_cell->x_min + 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min + 0.5*(curr_cell->y_max-curr_cell->y_min), curr_cell->y_max ); 
                        break;
                    case 1:
                        init_node(curr_cell->children[a], old_particle, curr_cell->x_min +0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                        break;
                    case 2:
                        init_node(curr_cell->children[a], old_particle, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                        break;
                    case 3:
                        init_node(curr_cell->children[a], old_particle, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min + 0.5*(curr_cell->y_max -curr_cell->y_min ), curr_cell->y_max ); 
                        break;
                    default: 
                        break;
                }
                //now old particle sits in correct node
            } else { //the current node has no particle. Then we're ok: 
                curr_cell->x_cm = (curr_cell->x_cm*curr_cell->M + m*x) /(m+curr_cell->M);
                curr_cell->y_cm = (curr_cell->y_cm*curr_cell->M + m*y) /(m+curr_cell->M);
                curr_cell->M = curr_cell->M + m;
            }
            //select the correct quadrant for the current particle
            b = getQuadrant(*it, curr_cell);
            //and adjust the boundaries
            adjust(b, x_min, x_max, y_min, y_max);

            //see if the quadrant that curr_particle should occupy is already there (meaning there are particles inside)
            if(!curr_cell->children[b]){ //if not, then create it and break the while loop. 
                init_node(curr_cell->children[b], &(*it), x_min, x_max, y_min, y_max);
                break;
            } else { //otherwise, just keep going
                curr_cell = curr_cell->children[b];
            }
        }
    }
    return ROOT;
}

void freeQuadTree(node* ROOT){
    if (ROOT == NULL) {
        return; 
    }
    for (int i = 0; i < 4; i++) {
        if (ROOT->children[i] != NULL) {
            freeQuadTree(ROOT->children[i]); 
        }
    }
    delete ROOT;
}

double width(node const*const n){
    return n->x_max-n->x_min;
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


//given a particle part, traverse the quadtree with root node n and compute the applied force usinf the BH approx
void calculate_force(Particle& part, node*& n, double theta){
    //std::cerr << "Particle: " << part.m;
    double forc[2];
    if (n->hasParticles) { //an external node
        if(n->p->id != part.id){ //a different particle
            force(forc, part, n->M, n->x_cm, n->y_cm);
            part.a_x += forc[0]/part.m; 
            part.a_y += forc[1]/part.m; 
        } else { //external node with the same particle
            part.a_x += 0; 
            part.a_y += 0; 
        }
    } else { //so an internal node
        if (width(n)/distance(part, n->x_cm, n->y_cm) < theta){ //
            force(forc, part, n->M, n->x_cm, n->y_cm);
            part.a_x += forc[0]/part.m; 
            part.a_y += forc[1]/part.m; 
        } else { 
            for (int i = 0; i < 4; i++){
                if(n->children[i]){
                    calculate_force(part, n->children[i], theta); 
                }
            }
        }
    }
}

void calculate_force(Particle& p1, Particle& p2){
    double forc[2];
    force(forc, p1, p2.m, p2.x, p2.y);
    p1.a_x += forc[0]/p1.m; 
    p1.a_y += forc[1]/p1.m; 
    p2.a_x += -forc[0]/p2.m; 
    p2.a_y += -forc[1]/p2.m; 
}


void leap_frog_step(Particle& part, double dt){
        part.x = part.x + part.v_x *dt; 
        part.y = part.y + part.v_y *dt; 
        part.v_x = part.v_x + part.a_x * dt; 
        part.v_y = part.v_y + part.a_y * dt; 
}

void barnesHuttAlgorithm(std::vector<Particle> & particles, double theta, double dt, int t){
    //first iteration using leap-frog
    if(t == 0){
        node* ROOT = nullptr;
        //build the quadtree
        ROOT = buildQuadTreeParallel(particles);
        //first step of leapfrog algorithm
        for (Particle& p : particles){
            calculate_force(p, ROOT, theta); 
            p.v_x = p.v_x + p.a_x * dt/2.; 
            p.v_y = p.v_y + p.a_y * dt/2.; 
            p.x = p.x + p.v_x * dt/2.; 
            p.y = p.y + p.v_y * dt/2.; 
        }
        freeQuadTree(ROOT);
    }

    // at every iteration, construct the quadtree
    node* ROOT = buildQuadTreeParallel(particles);
    // calculate forces using BH algo
    #pragma omp parallel for
        for (int i = 0; i < particles.size(); i++) {
            calculate_force(particles[i], ROOT, theta); 
            leap_frog_step(particles[i], dt);
        }
    //free the memory (to avoid memory leaks)
    freeQuadTree(ROOT);
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
        leap_frog_step(*it, dt);
    }
}