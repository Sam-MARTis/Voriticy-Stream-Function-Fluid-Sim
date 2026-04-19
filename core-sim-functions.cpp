#include "aux.hpp"
#include "core-sim-functions.hpp"
#include <algorithm> 
#include <cmath>
#include <iostream>

void advect_vorticity(float* ω, const float* x, const float* u, float u0, int nx, int ny, float dt, const float* dims);
void apply_viscosity(float* ω, int nx, int ny, float nu, float dt, const float* dims);

void solve_vorticity_transport(float* ω, const float* x, const float* u, const float u0, const float* ψ, const int nx, const int ny,const float nu, const float dt, const float* dims, const bool enable_advect_vorticity, const bool enable_apply_viscosity) {
    if(enable_advect_vorticity) {
        advect_vorticity(ω, x, u, u0, nx, ny, dt, dims);
    }
    if(enable_apply_viscosity) {
        apply_viscosity(ω, nx, ny, nu, dt, dims);
    }
    solve_boundary_vorticity_values(ω, u0, ψ, nx, ny, dims);
}

void advect_vorticity(float* ω, const float* x, const float* u, const float u0, const int nx, const int ny, const float dt, const float* dims) {
    float* ω_new = new float[nx * ny];
    std::cout<<"Advecting vorticity with dt: "<<dt<<std::endl;
     const float a = dims[0];
    for(int i=0; i < nx; i++) {
        for(int j=0; j < ny; j++) {
            const int idx = j*nx + i;
            const float px = x[2*idx];
            const float py = x[2*idx + 1];
            float u_x1, u_y1;
            find_velocity_at_point(u_x1, u_y1, px, py, u, u0, nx, ny, dims);
            float k1_x = -u_x1;
            float k1_y = -u_y1;
            float mid_px = px + 0.5f * dt * k1_x;
            float mid_py = py + 0.5f * dt * k1_y;
            float u_x2, u_y2;
            find_velocity_at_point(u_x2, u_y2, mid_px, mid_py, u, u0, nx, ny, dims);
            float k2_x = -u_x2;
            float k2_y = -u_y2;
            mid_px = px + 0.5f * dt * k2_x;
            mid_py = py + 0.5f * dt * k2_y; 
            float u_x3, u_y3;
            find_velocity_at_point(u_x3, u_y3, mid_px, mid_py, u, u0, nx, ny, dims);
            float k3_x = -u_x3;
            float k3_y = -u_y3;
            float end_px = px + dt * k3_x;
            float end_py = py + dt * k3_y;
            float u_x4, u_y4;
            find_velocity_at_point(u_x4, u_y4, end_px, end_py, u, u0, nx, ny, dims);
            float k4_x = -  u_x4;
            float k4_y = -  u_y4;
            float back_px = px + (dt / 6.0f) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
            float back_py = py + (dt / 6.0f) * (k1_y + 2*k2_y + 2*k3_y + k4_y);

            float ω_back;
            find_vorticity_at_point(ω_back, back_px, back_py, ω, nx, ny, dims);
            ω_new[idx] = ω_back;
        }
    }
    std::copy(ω_new, ω_new + (nx * ny), ω);
    delete[] ω_new;

}
void apply_viscosity(float* ω, const int nx, const int ny, const float nu, const float dt, const float* dims) {

    // const float inv_dx_squared = 1.0f / (dx * dx);
    // const float inv_dy_squared = 1.0f / (dy * dy);
    const float dξ = 1.0f/nx;
    const float dη = 1.0f/ny;
    const float inv_dξ_squared = 1.0f / (dξ * dξ);
    const float inv_dη_squared = 1.0f / (dη * dη);
    const float inv_dξdη = 1.0f / (dξ * dη);
    if(dt > 0.5f * (1.0f/(nu * (inv_dξ_squared + inv_dη_squared)))) {
        print("Warning: Time step may be too large for stability with the given viscosity.");
    }
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float ξx = 1.0f/a;
    const float ξy = 1.0f/(std::tanf(θ) * a);
    const float ηx = 0.0f;
    const float ηy = 1.0f/(std::sinf(θ) * b);

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



void solve_stream_function_update(float* ψ, const float* ω, const int nx, const int ny, const float* dims, const int max_iterations, const float tolerance) {
    // Simple Jacobi iteration for solving Poisson equation ∇²ψ = -ω
    float* ψ_new = new float[nx * ny];
    float* ψ_curr = ψ;
    float* ψ_next = ψ_new;
    for(int i=0; i<nx; i++){
        ψ_new[i] = 0.0f; // top boundary
        ψ[i] = 0.0f; 
        ψ_new[(ny-1)*nx + i] = 0.0f; // bottom boundary
        ψ[(ny-1)*nx + i] = 0.0f;
    }
    for(int j=0; j<ny; j++) {
        ψ_new[j*nx] = 0.0f; // left boundary
        ψ[j*nx] = 0.0f;
        ψ_new[j*nx + (nx-1)] = 0.0f; // right boundary
        ψ[j*nx + (nx-1)] = 0.0f;
    }
    // const float inv_dx_squared = 1.0f / (dx * dx);
    // const float inv_dy_squared = 1.0f / (dy * dy);
    const float dξ = 1.0f/nx;
    const float dη = 1.0f/ny;
    const float inv_dξ_squared = 1.0f / (dξ * dξ);
    const float inv_dη_squared = 1.0f / (dη * dη);
    const float inv_dξdη = 1.0f / (dξ * dη);
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float ξx = 1.0f/a;
    const float ξy = 1.0f/(std::tanf(θ) * a);
    const float ηx = 0.0f;
    const float ηy = 1.0f/(std::sinf(θ) * b);
    const float β = dξ/dη;
    const float a1 = ξx*ξx + ξy*ξy;
    const float a2 = (ξx*ηx + ξy*ηy) * 0.5f * β;
    const float a3 = (ηx*ηx + ηy*ηy) * β * β;

    for(int iter = 0; iter < max_iterations; iter++) {
        float max_diff = 0.0f;
        for(int i=1; i < nx-1; i++) {
            for(int j=1; j < ny-1; j++) {
                const int idx = j*nx + i;
                const float& ψ_left =  ψ_curr[j*nx + (i-1)];
                const float& ψ_right = ψ_curr[j*nx + (i+1)];
                const float& ψ_down = ψ_curr[(j-1)*nx + i];
                const float& ψ_up =  ψ_curr[(j+1)*nx + i];
                const float& ψ_down_left = ψ_curr[(j-1)*nx + (i-1)];
                const float& ψ_down_right = ψ_curr[(j-1)*nx + (i+1)];
                const float& ψ_up_left = ψ_curr[(j+1)*nx + (i-1)];
                const float& ψ_up_right = ψ_curr[(j+1)*nx + (i+1)];  
                const float rhs = -ω[idx];

                ψ_next[idx] = (1/(2*(a1 + a3)))*(a1*(ψ_left + ψ_right) + a2*(ψ_up_right - ψ_up_left - ψ_down_right + ψ_down_left) + a3*(ψ_down + ψ_up) - rhs*dξ*dξ);
                max_diff = std::max(max_diff, std::abs(ψ_next[idx] - ψ_curr[idx]));
            }
        }
        std::swap(ψ_curr, ψ_next);
        if(max_diff < tolerance) {
            break;
        }
    }

    if(ψ_curr != ψ) {
        std::copy(ψ_curr, ψ_curr + (nx * ny), ψ);
    }

    delete[] ψ_new;
}

void solve_velocity_update(float* u, const float u0, const float* ψ, const int nx, const int ny, const float* dims) {
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float ξx = 1.0f/a;
    const float ξy = 1.0f/(std::tanf(θ) * a);
    const float ηx = 0.0f;
    const float ηy = 1.0f/(std::sinf(θ) * b);
    const float dξ = 1.0f/nx;
    const float dη = 1.0f/ny;
    const float t1 = ξy/(2*dξ);
    const float t2 = ηy/(2*dη);
    const float t3 = ξx/(2*dξ);
    const float t4 = ηx/(2*dη);
    for(int i=1; i < nx-1; i++) {
        for(int j=1; j < ny-1; j++) {
            const int idx = j*nx + i;
            const float& ψ_left =  ψ[j*nx + (i-1)];
            const float& ψ_right = ψ[j*nx + (i+1)];
            const float& ψ_down = ψ[(j-1)*nx + i];
            const float& ψ_up =  ψ[(j+1)*nx + i];
            u[2*idx] = (ψ_right - ψ_left) * t1 + (ψ_up - ψ_down) * t2;
            u[2*idx + 1] = -((ψ_right - ψ_left) * t3 + (ψ_up - ψ_down) * t4);
        }
    }
    for(int i=0; i < nx; i++) {
        u[nx*(ny-1)*2 + 2*i] = u0; // top boundary
        u[nx*(ny-1)*2 + 2*i + 1] = 0.0f;
        u[i*2] = 0.0f; // bottom boundary
        u[i*2 + 1] = 0.0f;
    }
    for(int j=0; j < ny; j++) {
        u[j*nx*2] = 0.0f; // left boundary
        u[j*nx*2 + 1] = 0.0f;
        u[j*nx*2 + (nx-1)*2] = 0.0f; // right boundary
        u[j*nx*2 + (nx-1)*2 + 1] = 0.0f;
    }
}

void solve_boundary_vorticity_values(float* ω, const float u0, const float* ψ, const int nx, const int ny, const float* dims) {
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float ξx = 1.0f/a;
    const float ξy = 1.0f/(std::tanf(θ) * a);
    const float ηx = 0.0f;
    const float ηy = 1.0f/(std::sinf(θ) * b);
    const float dξ = 1.0f/nx;
    const float dη = 1.0f/ny;
    const float a1 = (ξx*ξx + ξy*ξy)/(dξ*dξ);
    const float a2 = (ξx*ηx + ξy*ηy)/(4*dξ*dη);
    const float a3 = (ηx*ηx + ηy*ηy)/(dη*dη);
    const float yη = std::sinf(θ) * b;
    for(int i=0; i < nx; i++) {
        const bool is_valid_horiz_walls  = (ψ[i] == 0.0f) && (ψ[(ny-1)*nx + i] == 0.0f);
        if(!is_valid_horiz_walls) {
            print("Error: Stream function values at horizontal boundaries must be zero for correct vorticity boundary conditions.");
            return;
        }
    }
    for(int j=0; j < ny; j++) {
        const bool is_valid_vert_walls  = (ψ[j*nx] == 0.0f) && (ψ[j*nx + (nx-1)] == 0.0f);
        if(!is_valid_vert_walls) {
            print("Error: Stream function values at vertical boundaries must be zero for correct vorticity boundary conditions.");
            return;
        }
    }

    for(int j=1; j < ny-1; j++) {
        // Wall 1
        ω[j*nx] = -2*(a1*ψ[j*nx + 1] + a2*(ψ[(j+1)*nx + 1] - ψ[(j-1)*nx + 1]));

        // Wall 3
        // const float d2ψdξdη_right = -(ψ[(j+1)*nx + (nx-2)] - ψ[(j-1)*nx + (nx-2)])/(4*dξ*dη);
        ω[j*nx + (nx-1)] = -2*(a1*ψ[j*nx + (nx-2)]+ a2*(-(ψ[(j+1)*nx + (nx-2)] - ψ[(j-1)*nx + (nx-2)])));
    }

    for(int i=1; i < nx-1; i++) {
        // Wall 2
        // const float d2ψdξdη_bottom = (ψ[nx + (i+1)] - ψ[nx + (i-1)] );
        ω[i] = -2*((ψ[nx + (i+1)] - ψ[nx + (i-1)] ) * a2 + a3*ψ[nx + i]);
        
        // Wall 4
        // const float d2ψdξdη_top = -(ψ[(ny-1)*nx + (i+1)] - ψ[(ny-1)*nx + (i-1)])/(4*dξ*dη);
        ω[(ny-1)*nx + i] = -2*(-(ψ[(ny-1)*nx + (i+1)] - ψ[(ny-1)*nx + (i-1)] ) * a2 + a3*(ψ[(ny-1)*nx + i] + u0*yη));

    }

    // Corner points
    ω[0] = -2*(a2*(ψ[nx + 1])); // Bottom-left corner
    ω[(ny-1)*nx] = -2*(-a2*ψ[((ny-1)-1)*nx + 1]); // Top-left corner
    ω[(nx-1)] = -2*(a2*(-(ψ[nx + (nx-2)]))); // Bottom-right corner
    ω[ny*nx -1] = -2*( a2*(ψ[((ny-1)-1)*nx + (nx-2)])); // Top-right corner
}