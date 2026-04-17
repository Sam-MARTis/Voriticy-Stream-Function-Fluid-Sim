#include "aux.hpp"
#include <algorithm> 
#include <cmath>
void solve_vorticity_transport(float* ω, const float* x, const float* u, const float u0, const int nx, const int ny, const float dx, const float dy, const float nu, const float dt, const float* dims) {
    advect_vorticity(ω, x, u, u0, nx, ny, dx, dy, dt, dims);
    apply_viscosity(ω, nx, ny, nu, dx, dy, dt, dims);
}

void advect_vorticity(float* ω, const float* x, const float* u, const float u0, const int nx, const int ny, const float dx, const float dy, const float dt, const float* dims) {
    for(int i=0; i < nx; i++) {
        for(int j=0; j < ny; j++) {
            const int idx = j*nx + i;
            const float px = x[2*idx];
            const float py = x[2*idx + 1];
            float u_x1, u_y1;
            find_velocity_at_point(u_x1, u_y1, px, py, u, u0, nx, ny, dx, dy);
            float k1_x = -u_x1;
            float k1_y = -u_y1;
            float mid_px = px + 0.5f * dt * k1_x;
            float mid_py = py + 0.5f * dt * k1_y;
            float u_x2, u_y2;
            find_velocity_at_point(u_x2, u_y2, mid_px, mid_py, u, u0, nx, ny, dx, dy);
            float k2_x = -u_x2;
            float k2_y = -u_y2;
            mid_px = px + 0.5f * dt * k2_x;
            mid_py = py + 0.5f * dt * k2_y; 
            float u_x3, u_y3;
            find_velocity_at_point(u_x3, u_y3, mid_px, mid_py, u, u0, nx, ny, dx, dy);
            float k3_x = -u_x3;
            float k3_y = -u_y3;
            float end_px = px + dt * k3_x;
            float end_py = py + dt * k3_y;
            float u_x4, u_y4;
            find_velocity_at_point(u_x4, u_y4, end_px, end_py, u, u0, nx, ny, dx, dy);
            float k4_x = -  u_x4;
            float k4_y = -  u_y4;
            float back_px = px + (dt / 6.0f) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
            float back_py = py + (dt / 6.0f) * (k1_y + 2*k2_y + 2*k3_y + k4_y);

            float ω_back;
            find_vorticity_at_point(ω_back, back_px, back_py, ω, nx, ny, dx, dy);
            ω[idx] = ω_back;
        }
    }
}
void apply_viscosity(float* ω, const int nx, const int ny, const float nu, const float dx, const float dy, const float dt, const float* dims) {

    const float inv_dx_squared = 1.0f / (dx * dx);
    const float inv_dy_squared = 1.0f / (dy * dy);
    const float dξ = 1.0f/nx;
    const float dη = 1.0f/ny;
    const float inv_dξ_squared = 1.0f / (dξ * dξ);
    const float inv_dη_squared = 1.0f / (dη * dη);
    const float inv_dξdη = 1.0f / (dξ * dη);
    if(dt > 0.5f * (1/(nu * (inv_dx_squared + inv_dy_squared)))) {
        print("Warning: Time step may be too large for stability with the given viscosity.");
    }
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float ξx = 1/a;
    const float ξy = 1/(std::tanf(θ) * a);
    const float ηx = 0.0f;
    const float ηy = 1/(std::sinf(θ) * b);

    for(int i=1; i < nx-1; i++) {
        for(int j=1; j < ny-1; j++) {
            const int idx = j*nx + i;
            const float& ω_center = ω[idx];
            const float& ω_left =  ω[j*nx + (i-1)];
            const float& ω_right = ω[j*nx + (i+1)];
            const float& ω_down = ω[(j-1)*nx + i];
            const float& ω_up =  ω[(j+1)*nx + i];
            const float& ω_down_left = ω[(j-1)*nx + (i-1)];
            const float& ω_down_right = ω[(j-1)*nx + (i+1)];
            const float& ω_up_left = ω[(j+1)*nx + (i-1)];
            const float& ω_up_right = ω[(j+1)*nx + (i+1)];  
            const float d2dξ2 = (ω_right - 2*ω_center + ω_left) * inv_dξ_squared;
            const float d2dη2 = (ω_up - 2*ω_center + ω_down) * inv_dη_squared;
            const float d2dξdη = (ω_up_right - ω_up_left - ω_down_right + ω_down_left) * 0.25f * inv_dξdη;
            const float laplacian = d2dξ2*(ξx*ξx + ξy*ξy) + d2dη2*(ηx*ηx + ηy*ηy) + 2 * d2dξdη * (ξx*ηx + ξy*ηy);

            ω[idx] += nu * laplacian * dt;
        }
    }
}
