import numpy as np
from src.solver import wave_step


def _stable_dt(c, dx, dy, safety=0.9):
    # CFL condition for the 2D wave equation (explicit):
    # (c*dt/dx)^2 + (c*dt/dy)^2 <= 1  ->  dt <= 1/(c*sqrt(1/dx^2 + 1/dy^2))
    denom = c * np.sqrt((1.0 / (dx * dx)) + (1.0 / (dy * dy)))
    return safety / denom


def _make_grid(nx, ny):
    x = np.arange(nx) - nx / 2.0
    y = np.arange(ny) - ny / 2.0
    return np.meshgrid(x, y, indexing="ij")


def _initial_conditions(nx, ny, X, Y, c_outside, mode="gaussian", amplitude=1.0, sigma=8.0, plane_xc=None, plane_lam=12.0):

    if mode == "gaussian":
        u0 = amplitude * np.exp(-(X**2 + Y**2) / (2.0 * sigma**2))
        v0 = np.zeros_like(u0)
        return u0, v0

    if mode == "ring":
        r = np.sqrt(X**2 + Y**2)
        radius = 12.0
        thickness = max(1.0, sigma / 2.0)
        u0 = amplitude * np.exp(-((r - radius) ** 2) / (2.0 * thickness**2))
        v0 = np.zeros_like(u0)
        return u0, v0

    if mode == "impulse_velocity":
        u0 = np.zeros((nx, ny))
        v0 = np.zeros((nx, ny))
        v0[nx // 2, ny // 2] = amplitude
        return u0, v0

    if mode == "plane_wave":
        if plane_xc is None:
            plane_xc = -nx / 3.0
        delta = X - plane_xc
        k = 2*np.pi/plane_lam
        u0 = amplitude*np.exp((-delta**2/(2*sigma**2)))*np.sin(k*delta)
        v0 = -c_outside*amplitude*np.exp((-delta**2/(2*sigma**2)))*(-delta/sigma**2*np.sin(k*delta)+k*np.cos(k*delta))
        return u0, v0

    raise ValueError(f"Unknown initial condition mode: {mode!r}")

def _make_c_field(X, Y, lens_x0, lens_y0, c_inside, c_outside, lens_R):
    r2 = (X-lens_x0)**2 + (Y-lens_y0)**2
    inside = (r2 <= lens_R**2) 
    c_grid = np.where(inside, c_inside, c_outside)
    return c_grid

def _make_gamma_field(nx, ny, width=15, gamma_max=2.0, power=2):
    if width <= 0 or gamma_max <= 0:
        return np.zeros((nx, ny), dtype=float)

    i = np.arange(nx)
    j = np.arange(ny)
    d_x = np.minimum(i, nx-1-i)
    d_y = np.minimum(j, ny-1-j)
    d_xy = np.minimum(d_x[:, None], d_y)
    s = (width - d_xy) / float(width)
    s = np.clip(s, 0.0, 1.0)
    gamma = gamma_max * (s ** power)
    return gamma

def _apply_bc_inplace(u, bc="neumann"):

    if bc == "neumann":
        u[0, :] = u[1, :]
        u[-1, :] = u[-2, :]
        u[:, 0] = u[:, 1]
        u[:, -1] = u[:, -2]
        return

    if bc == "dirichlet":
        u[0, :] = 0
        u[-1, :] = 0
        u[:, 0] = 0
        u[:, -1] = 0
        return

    if bc == "periodic":
        top = u[-2,:].copy()
        bottom = u[1, :].copy() 
        left = u[:, -2].copy()
        right = u[:, 1].copy()
        u[0, :] = top
        u[-1, :] = bottom
        u[:, 0] = left
        u[:, -1] = right
        return
    
    raise ValueError(f"Unknown bc: {bc!r}. Use 'neumann', 'dirichlet', or 'periodic'.")

def run_simulation(
    nx=120,
    ny=120,
    steps=3000,
    dx=1.0,
    dy=1.0,
    dt=None,
    record_every=5,
    init_mode="gaussian",
    init_amplitude=1.0,
    init_sigma=10.0,
    init_plane_xc=None,
    init_plane_lam=12.0,
    lens_x0=0,
    lens_y0=0,
    c_inside=6,
    c_outside=8,
    lens_R=20,
    bc="neumann",
    return_lens_mask=False,
    width=0,
    gamma_max=0.0,
    power=2,
):
    """Run a 2D wave simulation and return a list of snapshots (history)."""

    X, Y = _make_grid(nx, ny)    
    c_grid = _make_c_field(X, Y, lens_x0, lens_y0, c_inside, c_outside, lens_R)
    lens_mask = ((X - lens_x0)**2 + (Y - lens_y0)**2 <= lens_R**2)
    gamma = _make_gamma_field(nx, ny, width, gamma_max, power)

    if dx <= 0 or dy <= 0:
        raise ValueError("dx and dy must be > 0")

    if dt is None:
        dt = _stable_dt((c_grid.max()), dx, dy, safety=0.9)
    if dt <= 0:
        raise ValueError("dt must be > 0")
    if record_every <= 0:
        raise ValueError("record_every must be > 0")

    u, v0 = _initial_conditions(
        nx,
        ny,
        X,
        Y,
        c_outside,
        mode=init_mode,
        amplitude=init_amplitude,
        sigma=init_sigma,
        plane_xc=init_plane_xc,
        plane_lam=init_plane_lam,
    )

    # Encode initial velocity using: u_prev = u(t0 - dt) = u(t0) - dt*v(t0)
    u_prev = u - dt * v0

    history = []

    for step in range(steps):
        _apply_bc_inplace(u, bc)
        _apply_bc_inplace(u_prev, bc)
        u_new = wave_step(u, u_prev, c_grid, dx, dy, dt, gamma)
        u_prev = u
        _apply_bc_inplace(u_new, bc)
        u = u_new

        if step % record_every == 0:
            history.append(u.copy())

    if return_lens_mask: return history, lens_mask
    return history
