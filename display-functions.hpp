#pragma once 
#include<SFML/Graphics.hpp>

void render_velocities(const float* x, const float* u, const float normalization_factor, const float* origin, const float scaling[2], sf::RenderWindow& window);