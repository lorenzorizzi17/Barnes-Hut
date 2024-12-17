#include<SFML/Graphics.hpp>


void handleEvents(sf::RenderWindow& window){
    sf::Event event;
    while (window.pollEvent(event)) {
        if (event.type == sf::Event::Closed)
            window.close();
    }
    window.clear();
}


sf::Text loadFont(sf::Font& font){
    if (!font.loadFromFile("Ubuntu-Regular.ttf")) {
        throw std::runtime_error{"Font not found."};
    }
    sf::Text fpsText;
    fpsText.setFont(font);
    fpsText.setCharacterSize(15);  // Imposta la dimensione del testo
    fpsText.setFillColor(sf::Color::White);  // Imposta il colore del testo
    return fpsText;
}


