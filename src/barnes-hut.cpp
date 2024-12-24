#include"barnes-hut.hpp"


// the full Barnes-Hut algorithm
void barnesHuttAlgorithm(std::vector<Particle>& particles, double theta, double dt, int t){
    //first iteration using leap-frog
    if(t == 1){
        node* ROOT = nullptr;
        //build the quadtree
        ROOT = buildQuadTree(particles);
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
    node* ROOT = buildQuadTree(particles);
    // calculate forces using BH algo
    #pragma omp parallel for
        for (int i = 0; i < particles.size(); i++) {
            calculate_force(particles[i], ROOT, theta); 
            eulerStep(particles[i], dt);
            //maintenance
            particles[i].a_x = 0;
            particles[i].a_y = 0;
        }
    //free the memory (to avoid memory leaks)
    freeQuadTree(ROOT);
}

//initialize a node with a given particle
void initNode(node* & n , Particle* P, double xmin, double xmax, double ymin, double ymax){
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

// build the q-tree using the particle vector as input. Returns a pointer to the first element of the q-tree
node* buildQuadTree(std::vector<Particle>& particles){
    node* ROOT = NULL;
    //initialize the first cell with the fisrt particle
    initNode(ROOT, &particles.front(), -12, 12, -12, 12);

    //now run the algorithm for every particles
    auto it = particles.begin(); it++;
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

        // enter the loop until a break statement is found. Basically, we're traversing the q-tree
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
                        initNode(curr_cell->children[a], old_particle, curr_cell->x_min + 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min + 0.5*(curr_cell->y_max-curr_cell->y_min), curr_cell->y_max ); 
                        break;
                    case 1:
                        initNode(curr_cell->children[a], old_particle, curr_cell->x_min +0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                        break;
                    case 2:
                        initNode(curr_cell->children[a], old_particle, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                        break;
                    case 3:
                        initNode(curr_cell->children[a], old_particle, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min + 0.5*(curr_cell->y_max -curr_cell->y_min ), curr_cell->y_max ); 
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
                initNode(curr_cell->children[b], &(*it), x_min, x_max, y_min, y_max);
                break;
            } else { //otherwise, just keep going
                curr_cell = curr_cell->children[b];
            }
        }
    }
    return ROOT;
}

//this function splits the bodies into the four children and launch 4 tasks for every child. Used when building the q-tree parallel
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
                    initNode(curr_cell->children[a], p, curr_cell->x_min + 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min + 0.5*(curr_cell->y_max-curr_cell->y_min), curr_cell->y_max ); 
                    break;
                case 1:
                    initNode(curr_cell->children[a], p, curr_cell->x_min +0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->x_max,curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                    break;
                case 2:
                    initNode(curr_cell->children[a], p, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min, curr_cell->y_max - 0.5*(curr_cell->y_max-curr_cell->y_min) ); 
                    break;
                case 3:
                    initNode(curr_cell->children[a], p, curr_cell->x_min, curr_cell->x_max - 0.5*(curr_cell->x_max-curr_cell->x_min),curr_cell->y_min + 0.5*(curr_cell->y_max -curr_cell->y_min ), curr_cell->y_max ); 
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

// build the q-tree in a parallel fashion
node* buildQuadTreeParallel(std::vector<Particle>& particles){

    node* ROOT = NULL;
    initNode(ROOT, &*particles.begin(), -12, 12, -12, 12);
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

// starting from the root node, deletes the allocated memory of the q-tree
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


// debug function
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

// adjust the corner of the children according to parameter a
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

double width(node const*const n){
    return n->x_max-n->x_min;
}


//given a particle part, traverse the quadtree with root node n and compute the applied force usinf the BH approx
void calculate_force(Particle& part, node*& n, double theta){
    double forc[2];
    if (n->hasParticles) { //an external node
        if(n->p->id != part.id){ //a different particle
            force(forc, part, n->M, n->x_cm, n->y_cm);
            part.a_x += forc[0]/part.m; 
            part.a_y += forc[1]/part.m; 
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
