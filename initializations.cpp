#include<cmath>
void setup_inital_state(float* ψ, float* ω, float* x, float* u, const int nx, const int ny, const float dx, const float dy, const float θ, const float u0) {
    for (int i = 0; i < nx * ny; i++) {
        ψ[i] = 0.0f;
        ω[i] = 0.0f;
        u[2*i] = 0.0f;
        u[2*i + 1] = 0.0f;
        const float ξ = (i % nx) * dx;
        const float η = (i / nx) * dy;
        x[2*i] = ξ + η * std::cos(θ);
        x[2*i + 1] = η * std::sin(θ);
        u[2*i] = 0.0f;
        u[2*i + 1] = 0.0f;
    }
    for(int i = 0; i < nx; i++) {
        u[nx*(ny-1)*2 + 2*i] = u0;
        ω[nx*(ny-1) + i] = -2*u0 / dy;
    }
}

