def wave_step(u, u_prev, c_grid, dx, dy, dt, gamma):
    """Advance the 2D wave equation by one time step (explicit finite differences)."""

    c2_grid_x_left = (c_grid[2:, 1:-1]**2 + c_grid[1:-1, 1:-1]**2) / 2
    c2_grid_x_right = (c_grid[1:-1:, 1:-1]**2 + c_grid[:-2, 1:-1]**2) / 2
    c2_grid_y_left = (c_grid[1:-1, 2:]**2 + c_grid[1:-1, 1:-1]**2) / 2
    c2_grid_y_right = (c_grid[1:-1, 1:-1]**2 + c_grid[1:-1, :-2]**2) / 2

    F_x_left = c2_grid_x_left * (u[2:, 1:-1] - u[1:-1, 1:-1]) / dx
    F_x_right = c2_grid_x_right * (u[1:-1:, 1:-1] - u[:-2, 1:-1]) / dx
    F_y_left = c2_grid_y_left * (u[1:-1, 2:] - u[1:-1, 1:-1]) / dy
    F_y_right = c2_grid_y_right * (u[1:-1, 1:-1] - u[1:-1, :-2]) / dy
    
    div = (F_x_left - F_x_right) / dx + (F_y_left - F_y_right) / dy

    u_new = u.copy()
    u_new[1:-1, 1:-1] = 2 * u[1:-1, 1:-1] - u_prev[1:-1, 1:-1] + dt**2 * div - gamma[1:-1, 1:-1]*dt*(u[1:-1, 1:-1] - u_prev[1:-1, 1:-1])

    return u_new
