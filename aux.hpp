void print(const auto& message) {
    std::cout << message << std::endl;
}

void find_velocity_at_point(float& u_x, float& u_y, const float px, const float py, const float* u, const int nx, const int ny, const float dx, const float dy) {
    const int i = static_cast<int>(px / dx);
    const int j = static_cast<int>(py / dy);
    const float i_frac = (px - i * dx) / dx;
    const float j_frac = (py - j * dy) / dy;
    float ux[4] = {0, 0, 0, 0}; // bottom left, top left, top right, bottom right
    float uy[4] = {0, 0, 0, 0};
    const bool bottom_left_in_bound = (i >= 0 && i < nx-1) && (j >= 0 && j < ny-1);
    const bool bottom_right_in_bound = (i >= -1 && i < nx-2) && (j >= 0 && j < ny-1);
    const bool top_left_in_bound = (i >= 0 && i < nx-1) && (j >= -1 && j < ny-2);
    const bool top_right_in_bound = (i >= -1 && i < nx-2) && (j >= -1 && j < ny-2);
    const bool bounds[4] = {bottom_left_in_bound, top_left_in_bound, top_right_in_bound, bottom_right_in_bound};
    if (bounds[0]) {
        ux[0] = u[2*(j*nx + i)];
        uy[0] = u[2*(j*nx + i) + 1];
    }
    if (bounds[1]) {
        ux[1] = u[2*((j+1)*nx + i)];
        uy[1] = u[2*((j+1)*nx + i) + 1];
    }
    if (bounds[2]) {
        ux[2] = u[2*((j+1)*nx + (i+1))];
        uy[2] = u[2*((j+1)*nx + (i+1)) + 1];
    }
    if (bounds[3]) {
        ux[3] = u[2*(j*nx + (i+1))];
        uy[3] = u[2*(j*nx + (i+1)) + 1];
    }
    u_x = (ux[0] * (1-i_frac) + ux[3] * i_frac) * (1-j_frac) + (ux[1] * (1-i_frac) + ux[2] * i_frac) * j_frac;
    u_y = (uy[0] * (1-i_frac) + uy[3] * i_frac) * (1-j_frac) + (uy[1] * (1-i_frac) + uy[2] * i_frac) * j_frac;
}