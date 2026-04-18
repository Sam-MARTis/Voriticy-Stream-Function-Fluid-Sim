#include<iostream>
#include<cmath>
#include<SFML/Graphics.hpp>
#include "constants.hpp"
#include "core-sim-functions.hpp"
#include "display-functions.hpp"
#include "initializations.hpp"




void main() {
    const float a = 1.0f;
    const float b = 1.0f;
    const float θ = M_PI/2;
    const int nx = 100;
    const int ny = 100;
    // const float dx = ;
    // const float dy = 0.01f;

    const float nu = 0.1f;
    const float dt = 0.01f;
    const float u0 = 1.0f;
    const float dims[3] = {1.0f, 1.0f, M_PI/4};
    float* ψ = new float[nx*ny];
    float* ω = new float[nx*ny];
    float* u = new float[2*nx*ny];
    float* x = new float[2*nx*ny];





    setup_inital_state(ψ, ω, x, u, nx, ny, dims, u0);
    
    for(int iter=0; iter < 100; iter++) {
        solve_vorticity_transport(ω, x, u, u0, ψ, nx, ny, nu, dt, dims);
        solve_stream_function_update(ψ, ω, nx, ny, dims, 10000, 1e-6f);
        solve_velocity_update(u, u0, ψ, nx, ny, dims);
        solve_boundary_vorticity_values(ω, u0, ψ, nx, ny, dims);
        render_velocities(x, u, u0, dims,);
        if(iter % 10 == 0) {
            std::cout << "Iteration: " << iter << std::endl;
        }
    }
    delete[] ψ;
    delete[] ω;
    delete[] u;
}