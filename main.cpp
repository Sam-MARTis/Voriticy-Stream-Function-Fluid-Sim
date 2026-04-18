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
const int NX = 100;
const int NY = 100;
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
    const float dt = 0.001f;
    const float u0 = 1.0f;
    const float dims[3] = {a, b, θ};
    float* ψ = new float[NX*NY];
    float* ω = new float[NX*NY];
    float* u = new float[2*NX*NY];
    float* x = new float[2*NX*NY];


    window.create(sf::VideoMode({SCREEN_WIDTH + (SCREEN_OFFSET_X + SCREEN_END_X_PADDING), SCREEN_HEIGHT + (SCREEN_OFFSET_Y + SCREEN_END_Y_PADDING)}, 10), "Fluid Simulation");
    window.setFramerateLimit(FRAME_RATE_LIMIT);
    ImGui::SFML::Init(window);
    const float scaling[2] = {SCREEN_WIDTH/(a + b*std::cos(θ))*0.5, SCREEN_WIDTH/(a + b*std::cos(θ))*0.5};



    setup_inital_state(ψ, ω, x, u, NX, NY, dims, u0);
    
    for(int iter=0; iter < 100; iter++) {
        solve_vorticity_transport(ω, x, u, u0, ψ, NX, NY, nu, dt, dims);
        solve_stream_function_update(ψ, ω, NX, NY, dims, 10000, 1e-6f);
        solve_velocity_update(u, u0, ψ, NX, NY, dims);
        solve_boundary_vorticity_values(ω, u0, ψ, NX, NY, dims);
        render_velocities(x, u, u0, dims, scaling, window);
        if(iter % 10 == 0) {
            std::cout << "Iteration: " << iter << std::endl;
        }
    }
    delete[] ψ;
    delete[] ω;
    delete[] u;
    return 0;
}