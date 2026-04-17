#include "aux.hpp"

void solve_vorticity_transport(float* ω, const float* x, const float* u, const int nx, const int ny, const float dx, const float dy, const float nu, const float dt) {
    advect_vorticity(ω, x, u, nx, ny, dx, dy, dt);
}
void advect_vorticity(float* ω, const float* x, const float* u, const int nx, const int ny, const float dx, const float dy, const float dt) {
}
void apply_viscosity(float* ω, const int nx, const int ny, const float nu, const float dx, const float dy, const float dt) {
}