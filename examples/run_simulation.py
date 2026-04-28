import sys
import os
import argparse

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.simulation import run_simulation
from src.visualization import animate_simulation

parser = argparse.ArgumentParser(description="2D wave equation simulation (finite differences).")
parser.add_argument("--nx", type=int, default=120)
parser.add_argument("--ny", type=int, default=120)
parser.add_argument("--steps", type=int, default=3000)
parser.add_argument("--c-inside", "--c_inside", dest="c_inside", type=float, default=6)
parser.add_argument("--c-outside", "--c_outside", dest="c_outside", type=float, default=8)
parser.add_argument("--lens-x0", "--lens_x0", dest="lens_x0", type=float, default=0)
parser.add_argument("--lens-y0", "--lens_y0", dest="lens_y0", type=float, default=0)
parser.add_argument("--lens-r", "--lens_R", "--lens_r", dest="lens_R", type=float, default=20)
parser.add_argument("--dt", type=float, default=None)
parser.add_argument("--record-every", "--record_every", dest="record_every", type=int, default=5)
parser.add_argument(
    "--bc",
    dest="bc",
    type=str,
    default="neumann",
    choices=["neumann", "dirichlet", "periodic"],
)
parser.add_argument(
    "--init-mode",
    "--init_mode",
    dest="init_mode",
    type=str,
    default="gaussian",
    choices=["gaussian", "ring", "impulse_velocity", "plane_wave"],
)
parser.add_argument("--amplitude", type=float, default=1.0)
parser.add_argument("--sigma", type=float, default=10.0)
parser.add_argument(
    "--plane-xc",
    "--plane_xc",
    dest="plane_xc",
    type=float,
    default=None,
    help="Plane-wave packet center x (grid coordinates, centered at 0). Only used with --init-mode plane_wave.",
)
parser.add_argument(
    "--plane-lam",
    "--plane_lam",
    dest="plane_lam",
    type=float,
    default=12.0,
    help="Plane-wave wavelength in grid cells. Only used with --init-mode plane_wave.",
)
parser.add_argument(
    "--sponge-width",
    "--sponge_width",
    dest="sponge_width",
    type=int,
    default=0,
    help="Absorbing boundary (sponge layer) thickness in grid cells. 0 disables damping.",
)
parser.add_argument(
    "--sponge-gamma-max",
    "--sponge_gamma_max",
    dest="sponge_gamma_max",
    type=float,
    default=0.0,
    help="Max damping strength in the sponge layer. 0 disables damping.",
)
parser.add_argument(
    "--sponge-power",
    "--sponge_power",
    dest="sponge_power",
    type=int,
    default=2,
    help="Sponge ramp exponent (2 or 3 are common). Only used when sponge is enabled.",
)
parser.add_argument(
    "--save",
    dest="save",
    type=str,
    default=None,
    help="Save the animation to a file (e.g. out.gif or out.mp4).",
)
parser.add_argument(
    "--no-show",
    "--no_show",
    dest="no_show",
    action="store_true",
    help="Do not open an interactive window (useful together with --save).",
)
parser.add_argument(
    "--fps",
    dest="fps",
    type=float,
    default=None,
    help="Frames per second when saving (defaults to 1000/interval_ms).",
)
parser.add_argument(
    "--dpi",
    dest="dpi",
    type=int,
    default=120,
    help="DPI used when saving the animation.",
)
parser.add_argument("--interval-ms", "--interval_ms", dest="interval_ms", type=int, default=30)
args = parser.parse_args()

history, lens_mask = run_simulation(
    nx=args.nx,
    ny=args.ny,
    steps=args.steps,
    dt=args.dt,
    c_inside=args.c_inside,
    c_outside=args.c_outside,
    lens_x0=args.lens_x0,
    lens_y0=args.lens_y0,
    lens_R=args.lens_R,
    record_every=args.record_every,
    bc=args.bc,
    init_mode=args.init_mode,
    init_amplitude=args.amplitude,
    init_sigma=args.sigma,
    init_plane_xc=args.plane_xc,
    init_plane_lam=args.plane_lam,
    width=args.sponge_width,
    gamma_max=args.sponge_gamma_max,
    power=args.sponge_power,
    return_lens_mask=True
)

animate_simulation(
    history,
    interval_ms=args.interval_ms,
    lens_mask=lens_mask,
    save_path=args.save,
    show=not args.no_show,
    fps=args.fps,
    dpi=args.dpi,
)
