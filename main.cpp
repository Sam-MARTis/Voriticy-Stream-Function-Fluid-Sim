#include<iostream>
#include<cmath>
#include<SFML/Graphics.hpp>
#include "constants.hpp"
#include "core-sim-functions.hpp"
#include "display-functions.hpp"
#include "initializations.hpp"
#include "imgui.h"
#include "imgui-SFML.h"

sf::RenderWindow window;
const int NX = 30;
const int NY = 30;
const int SCREEN_WIDTH = SCREEN_WIDTH_default;
const int SCREEN_HEIGHT = SCREEN_HEIGHT_default;
const int SCREEN_OFFSET_X = SCREEN_OFFSET_X_default;
const int SCREEN_OFFSET_Y = SCREEN_OFFSET_Y_default;
const int SCREEN_END_X_PADDING = SCREEN_END_X_PADDING_default;
const int SCREEN_END_Y_PADDING = SCREEN_END_Y_PADDING_default;

int main() {
    const float a = 10.0f;
    const float b = 10.0f;
    const float θ = M_PI/2;

    // const int NX = 100;
    // const int NY = 100;
    // const float dx = ;
    // const float dy = 0.01f;

    const float nu = 0.1f;
    const float dt = 0.0001f;
    const float u0 = 1.0f;
    const float origin[2] = {static_cast<float>(SCREEN_OFFSET_X), static_cast<float>(SCREEN_OFFSET_Y)};
    const float dims[3] = {a, b, θ};
    float* ψ = new float[NX*NY];
    float* ω = new float[NX*NY];
    float* u = new float[2*NX*NY];
    float* x = new float[2*NX*NY];


    window.create(sf::VideoMode({SCREEN_WIDTH + (SCREEN_OFFSET_X + SCREEN_END_X_PADDING), SCREEN_HEIGHT + (SCREEN_OFFSET_Y + SCREEN_END_Y_PADDING)}, 10), "Fluid Simulation");
    window.setFramerateLimit(FRAME_RATE_LIMIT);
    ImGui::SFML::Init(window);
    const float scaling[2] = {SCREEN_WIDTH/(a + b*std::cos(θ))*0.5, SCREEN_WIDTH/(a + b*std::cos(θ))*0.5};

    int render_mode = 0;
    bool render_velocities_enabled = true;
    float normalization_constant = 1.0f;
    int velocity_thickness = 2;
    float low_colour[3] = {0.0f, 0.0f, 1.0f};
    float high_colour[3] = {1.0f, 0.0f, 0.0f};
    int iter = 0;
    sf::Clock deltaClock;



    setup_inital_state(ψ, ω, x, u, NX, NY, dims, u0);
    
    while(window.isOpen()) {
        while(const auto event = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *event);
            if(event->is<sf::Event::Closed>()) {
                window.close();
            }
        }

        ImGui::SFML::Update(window, deltaClock.restart());

        ImGui::Begin("Render Controls");
        ImGui::Text("Iteration: %d", iter);
        const char* render_modes[] = {"None", "Vorticity", "Stream Function"};
        ImGui::Combo("Field", &render_mode, render_modes, IM_ARRAYSIZE(render_modes));
        ImGui::ColorEdit3("Low Colour", low_colour);
        ImGui::ColorEdit3("High Colour", high_colour);
        ImGui::Checkbox("Render Velocities", &render_velocities_enabled);
        ImGui::SliderFloat("Normalization", &normalization_constant, 0.01f, 20.0f, "%.3f");
        ImGui::SliderInt("Thickness", &velocity_thickness, 1, 10);
        ImGui::End();

        solve_vorticity_transport(ω, x, u, u0, ψ, NX, NY, nu, dt, dims);

        // solve_stream_function_update(ψ, ω, NX, NY, dims, 10000, 1e-6f);

        // solve_velocity_update(u, u0, ψ, NX, NY, dims);

        // solve_boundary_vorticity_values(ω, u0, ψ, NX, NY, dims);

        window.clear(sf::Color::Black);

        const sf::Color low_sf(static_cast<unsigned char>(low_colour[0] * 255.0f),
                               static_cast<unsigned char>(low_colour[1] * 255.0f),
                               static_cast<unsigned char>(low_colour[2] * 255.0f));
        const sf::Color high_sf(static_cast<unsigned char>(high_colour[0] * 255.0f),
                                static_cast<unsigned char>(high_colour[1] * 255.0f),
                                static_cast<unsigned char>(high_colour[2] * 255.0f));

        if(render_mode == 1) {
            render_scalar_field(x, ω, NX, NY, origin, scaling, low_sf, high_sf, window);
        } else if(render_mode == 2) {
            render_scalar_field(x, ψ, NX, NY, origin, scaling, low_sf, high_sf, window);
        }

        if(render_velocities_enabled) {
            render_velocities(x, u, NX, NY, normalization_constant, velocity_thickness, origin, scaling, window);
        }

        ImGui::SFML::Render(window);
        window.display();

        if(iter % 10 == 0) {
            std::cout << "Iteration: " << iter << std::endl;
        }
        std::cout << "u at center: " << u[static_cast<int>((NX/2)*NY*2)] << std::endl;
        iter++;
    }

    ImGui::SFML::Shutdown();
    delete[] ψ;
    delete[] ω;
    delete[] u;
    delete[] x;
    return 0;
}