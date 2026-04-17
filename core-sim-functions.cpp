#include "aux.hpp"

void solve_vorticity_transport(float* ω, const float* x, const float* u, const float u0, const int nx, const int ny, const float dx, const float dy, const float nu, const float dt) {
    advect_vorticity(ω, x, u, u0, nx, ny, dx, dy, dt);
    apply_viscosity(ω, nx, ny, nu, dx, dy, dt);
}

void advect_vorticity(float* ω, const float* x, const float* u, const float u0, const int nx, const int ny, const float dx, const float dy, const float dt) {
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
void apply_viscosity(float* ω, const int nx, const int ny, const float nu, const float dx, const float dy, const float dt) {
    const float inv_dx_squared = 1.0f / (dx * dx);
    const float inv_dy_squared = 1.0f / (dy * dy);
    if(dt > 0.5f * (1/(nu * (inv_dx_squared + inv_dy_squared)))) {
        print("Warning: Time step may be too large for stability with the given viscosity.");
    }
    for(int i=1; i < nx-1; i++) {
        for(int j=1; j < ny-1; j++) {
            const int idx = j*nx + i;
            float ω_center = ω[idx];
            float ω_left =  ω[j*nx + (i-1)];
            float ω_right = ω[j*nx + (i+1)];
            float ω_down = ω[(j-1)*nx + i];
            float ω_up =  ω[(j+1)*nx + i];
            float laplacian = (ω_left - 2*ω_center + ω_right) * inv_dx_squared + (ω_down - 2*ω_center + ω_up) * inv_dy_squared;
            ω[idx] += nu * laplacian * dt;
        }
    }
}
