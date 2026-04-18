#include<SFML/Graphics.hpp>
#include <cmath>

void render_velocities(const float* x, const float* u, const float normalization_factor, const float* origin, const float* scaling, sf::RenderWindow& window) {
    const float inv_normalization_factor = 1.0f / normalization_factor;
    for(int i=0; i < 100; i++) {
        float px = origin[0] + x[2*i] * scaling[0];
        float py = origin[1] + x[2*i + 1] * scaling[1];
        float u_x = u[2*i] * inv_normalization_factor;
        float u_y = u[2*i + 1] * inv_normalization_factor;
        draw_arrow(px, py, u_x, u_y, window);
     }
}
void draw_arrow(const float px, const float py, const float dPx, const float dPy, sf::RenderWindow& window, const int thickness=2, const float head_fraction = 0.2,  sf::Color colour = sf::Color::Red) {
    sf::Vector2f start(px, py);
    sf::Vector2f direction(dPx, dPy);
    sf::Vector2f end = start + direction;
    const float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
    const sf::Angle angleDegrees = std::atan2(direction.y, direction.x) * sf::radians(1);
    const float head_length = length * head_fraction;
    sf::RectangleShape line(sf::Vector2f(length - head_length, (float)thickness));
    line.setOrigin({0, (float)thickness * 0.5f});
    line.setPosition(start);
    line.setRotation(angleDegrees);
    line.setFillColor(colour);

    sf::ConvexShape head;
    head.setPointCount(3);
    head.setPoint(0, {0, 0});
    head.setPoint(1, {-head_length, head_length * 0.5f});
    head.setPoint(2, {-head_length, -head_length * 0.5f});
    head.setFillColor(colour);

    head.setPosition(end);
    head.setRotation(angleDegrees);

    window.draw(line);
    window.draw(head);
}

// void drawArrow(sf::RenderWindow &window, sf::Vector2f start, sf::Vector2f end,
//                int thickness = 2, float head_fraction = 0.2f,
//                sf::Color colour = sf::Color::Red)
// {
//     sf::Vector2f direction = end - start;
//     float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
//     sf::Angle angleDegrees = std::atan2(direction.y, direction.x) * sf::radians(1);

//     float head_length = length * head_fraction;
//     sf::RectangleShape line(sf::Vector2f(length - head_length, (float)thickness));
//     line.setOrigin({0, (float)thickness * 0.5f});
//     line.setPosition(start);
//     line.setRotation(angleDegrees);
//     line.setFillColor(colour);

//     sf::ConvexShape head;
//     head.setPointCount(3);
//     head.setPoint(0, {0, 0});
//     head.setPoint(1, {-head_length, head_length * 0.5f});
//     head.setPoint(2, {-head_length, -head_length * 0.5f});
//     head.setFillColor(colour);

//     head.setPosition(end);
//     head.setRotation(angleDegrees);

//     window.draw(line);
//     window.draw(head);
// }
