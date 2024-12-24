#include "graphics.hpp"
#include "initial_distr.hpp"
#include "barnes-hut.hpp"
#include "nbody-sim.hpp"
#include <sstream>
#include <iomanip>

void Simulation::init(initialDistribution id, int N){
    switch (id)
    {
    case SINGLEGALAXY:
        makeRealisticGalaxy(m_particles, N);
        break;
    case DOUBLEGALAXY:
        makeDoubleGalaxies(m_particles, N/2, N/2);
    default:
        break;
    }
}

void Simulation::run(simulationMode mode, double dt){

    omp_set_num_threads(NUM_THREADS);
    //graphical parameters
    int DIM_ON_SCREEN = 1000;
    int DIM = 12;

    //SFML stuff: create a window, load the font and instantiate a clock
    sf::RenderWindow window = sf::RenderWindow(sf::VideoMode(DIM_ON_SCREEN, DIM_ON_SCREEN), "Barnes-Hut");
    sf::Font font;
    sf::Text fpsText = loadFont(font);
    sf::Text nText = loadFont(font);
    sf::Clock clock;
    sf::Clock general_clock;

    //graphical window and main loop
    while (window.isOpen() && m_time < 10000) {

        // -- handle the events and clear the window --
        handleEvents(window);

        // -- draw every particle --
        std::for_each(m_particles.begin(),m_particles.end(),[&](Particle const& p){drawParticle(window, p, DIM, DIM_ON_SCREEN);});

        // -- dynamics --
        switch (mode) {
            case BARNESHUTT:
                barnesHuttAlgorithm(m_particles,0.5 ,dt , m_time);
                break;
            case DOUBLELOOP:
                nestedLoopAlgorithm(m_particles, dt, m_time);
                break;
            default:
                break;
        }
        
        // -- draw the fps and display the window --
        std::ostringstream str;
        str << std::fixed << std::setprecision(2) << double(1/(clock.restart().asSeconds()));
        fpsText.setString("FPS = " + str.str());
        nText.setString("N = " + std::to_string(m_particles.size()));
        nText.setPosition(window.getSize().x - nText.getLocalBounds().width-1,0);
        window.draw(fpsText);
        window.draw(nText);
        window.display();

        for(auto it = m_particles.begin(); it != m_particles.end();){
            if(it->x <= -DIM || it->x >= DIM || it->y >= DIM || it->y <= -DIM){
                it = m_particles.erase(it);
            } else {
                it++;
            }
        }
        m_time++;
    }
    
}