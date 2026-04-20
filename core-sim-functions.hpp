#pragma once
void solve_vorticity_transport(float* ω, const float* x, const float* u, const float u0, const float* ψ, const int nx, const int ny,const float nu, const float dt, const float* dims, const bool use_operator_splitting, const bool enable_apply_viscosity, const bool enable_advect_vorticity);
void solve_stream_function_update(float* ψ, const float* ω, const int nx, const int ny, const float* dims, const int max_iterations, const float tolerance);
bool check_stream_function_convergence(const float* ψ, const float* ω, const int nx, const int ny, const float* dims, const float tolerance, float& max_residual);
void solve_velocity_update(float* u, const float u0, const float* ψ, const int nx, const int ny, const float* dims);
void solve_boundary_vorticity_values(float* ω, const float u0, const float* ψ, const int nx, const int ny, const float* dims);