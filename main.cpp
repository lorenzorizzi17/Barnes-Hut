#include"barnes-hutt.hpp"
#include"graphics.hpp"
#include<sstream>
#include<iomanip>

#define NUM_THREADS 16


int main(){

    //main parameters
    double theta = 0.5;
    double dt = 0.005;
    const int N = 100000;
    omp_set_num_threads(NUM_THREADS);
    //graphical parameters
    int DIM_ON_SCREEN = 1000;
    int DIM = 12;

    //particles are collected in a std::vector 
    std::vector<Particle> particles;
    particles.reserve(N+1);

    //initial condition: create a circular-ish distribution
    makeRealisticGalaxy(particles, N);
    //makeDoubleGalaxies(particles, N/2,N/2 );

    //SFML stuff: create a window, load the font and instantiate a clock
    sf::RenderWindow window = sf::RenderWindow(sf::VideoMode(DIM_ON_SCREEN, DIM_ON_SCREEN), "Barnes-Hut");
    sf::Font font;
    sf::Text fpsText = loadFont(font);
    sf::Text nText = loadFont(font);
    sf::Clock clock;
    sf::Clock general_clock;

    int time = 0;
    //graphical window and main loop
    while (window.isOpen() && time < 10000) {

        // -- handle the events and clear the window --
        handleEvents(window);

        // -- draw every particle --
        for (const Particle& p : particles) {
            sf::CircleShape point(std::log10(10*p.m));
            point.setOrigin(0.5*std::log10(10*p.m), 0.5*std::log10(10*p.m));
            point.setFillColor(sf::Color(255,255,255,125));
            float x_screen = ((p.x + DIM) / (DIM + DIM)) * DIM_ON_SCREEN; 
            float y_screen = DIM_ON_SCREEN - ((p.y + DIM ) / (DIM + DIM)) * DIM_ON_SCREEN; 
            point.setPosition(x_screen, y_screen);
            window.draw(point);  
        }

        // -- dynamics --
        barnesHuttAlgorithm(particles, theta, dt, time);
        //nestedLoopAlgorithm(particles,dt,time );
        
        // -- draw the fps and display the window --
        std::ostringstream str;
        str << std::fixed << std::setprecision(2) << double(1/(clock.restart().asSeconds()));
        fpsText.setString("FPS = " + str.str());
        nText.setString("N = " + std::to_string(particles.size()));
        nText.setPosition(window.getSize().x - nText.getLocalBounds().width-1,0);
        window.draw(fpsText);
        window.draw(nText);
        window.display();


        // -- set acc to zero and remove particles that crossed the borders --
        for(auto it = particles.begin(); it != particles.end(); it++){
            it->a_x = 0;
            it->a_y = 0;
        }
        for(auto it = particles.begin(); it != particles.end();){
            if(it->x <= -DIM || it->x >= DIM || it->y >= DIM || it->y <= -DIM){
                it = particles.erase(it);
            } else {
                it++;
            }
        }
        time++;
    }
        
    std::cout << general_clock.getElapsedTime().asSeconds() << " s\n";
}

