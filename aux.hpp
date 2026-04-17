void print(const auto& message);
void find_velocity_at_point(float& u_x, float& u_y, const float px, const float py, const float* u, const float u0, const int nx, const int ny, const float dx, const float dy);
void find_vorticity_at_point(float& ω_val, const float px, const float py, const float* ω, const int nx, const int ny, const float dx, const float dy);
