#ifndef GRAPHICS
#define GRAPHICS

#include<SFML/Graphics.hpp>
#include"gravity.hpp"

void handleEvents(sf::RenderWindow& window){
    sf::Event event;
    while (window.pollEvent(event)) {
        if (event.type == sf::Event::Closed)
            window.close();
    }
    window.clear();
}


sf::Text loadFont(sf::Font& font){
    if (!font.loadFromFile("res/Ubuntu-Regular.ttf")) {
        throw std::runtime_error{"Font not found."};
    }
    sf::Text fpsText;
    fpsText.setFont(font);
    fpsText.setCharacterSize(15);
    fpsText.setFillColor(sf::Color::White);  
    return fpsText;
}

void drawParticle(sf::RenderWindow& window, Particle const& p, int DIM, int DIM_ON_SCREEN){
            sf::CircleShape point(std::log10(10*p.m));
            point.setOrigin(0.5*std::log10(10*p.m), 0.5*std::log10(10*p.m));
            point.setFillColor(sf::Color(255,255,255,125));
            float x_screen = ((p.x + DIM) / (DIM + DIM)) * DIM_ON_SCREEN; 
            float y_screen = DIM_ON_SCREEN - ((p.y + DIM ) / (DIM + DIM)) * DIM_ON_SCREEN; 
            point.setPosition(x_screen, y_screen);
            window.draw(point);  
}


#endif