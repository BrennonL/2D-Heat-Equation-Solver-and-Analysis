import numpy as np
import matplotlib.pyplot as plt


# Exact solution of the 2D heat equation
def analytical_solution(X, Y, t, alpha, Lx=1.0, Ly=1.0, T0=100):
    return T0 * np.exp(-alpha * np.pi**2 * ((1 / Lx**2) + (1 / Ly**2)) * t) * \
           np.sin(np.pi * X / Lx) * np.sin(np.pi * Y / Ly)


# Explicit FTCS solver for the 2D transient heat equation.
def ftcs_2d_heat(
    Lx=1.0,
    Ly=1.0,
    nx=41,
    ny=41,
    alpha=1.0,
    dt=1e-5,
    t_final=0.01,
    T0=100
):
    """
    PDE:
        dT/dt = alpha * (d2T/dx2 + d2T/dy2)

    Boundary condition:
        T = 0 on all boundaries

    Initial condition:
        T(x,y,0) = T0*sin(pi*x)*sin(pi*y)
    """

    # Grid setup
    dx = Lx / (nx - 1) 
    dy = Ly / (ny - 1)

    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    X, Y = np.meshgrid(x, y)

    # Stability parameters
    rx = alpha * dt / dx**2
    ry = alpha * dt / dy**2

    # Check FTCS stability condition
    if rx + ry > 0.5:
        raise ValueError(
            f"FTCS stability condition violated: rx + ry = {rx + ry:.4f} > 0.5"
        )

    # Initial condition
    T = T0 * np.sin(np.pi * X / Lx) * np.sin(np.pi * Y / Ly)

    # Zero Dirichlet boundary conditions (T = 0 on all edges)
    T[0, :] = 0      # bottom
    T[-1, :] = 0     # top
    T[:, 0] = 0      # left
    T[:, -1] = 0     # right

    # Number of time steps
    n_steps = int(t_final / dt)

    # Time stepping loop (FTCS)
    for _ in range(n_steps):

        # Store previous time step
        T_old = T.copy() 

        # Update interior nodes using FTCS formula
        T[1:-1, 1:-1] = (
            T_old[1:-1, 1:-1]
            # second derivative in x-direction
            + rx * (T_old[1:-1, 2:] - 2*T_old[1:-1, 1:-1] + T_old[1:-1, :-2])
            # second derivative in y-direction
            + ry * (T_old[2:, 1:-1] - 2*T_old[1:-1, 1:-1] + T_old[:-2, 1:-1])
        )

        # Reapply boundary conditions
        T[0, :] = 0
        T[-1, :] = 0
        T[:, 0] = 0
        T[:, -1] = 0

    return x, y, X, Y, T, dx, dy, rx, ry


# Compute L2 and maximum error between numerical and exact solutions
def compute_errors(T_num, T_exact):
    error = T_num - T_exact

    L2_error = np.sqrt(np.mean(error**2))
    max_error = np.max(np.abs(error))

    return L2_error, max_error, error

# Plot temperature distribution as a heatmap
def plot_heatmap(X, Y, T, title):
    plt.figure(figsize=(6, 5))
    plt.contourf(X, Y, T, levels=30, cmap="hot")
    plt.colorbar(label="Temperature")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.tight_layout()


# Print numerical vs exact values at selected points
def print_selected_values(x, y, T_num, T_exact):
    points = [
        (0.25, 0.25),
        (0.50, 0.50),
        (0.75, 0.75),
        (0.25, 0.75),
        (0.75, 0.25)
    ]

    print("\nSelected Grid Point Comparison")
    print("-" * 75)
    print(f"{'x':>8} {'y':>8} {'FTCS':>15} {'Exact':>15} {'Abs Error':>15}")
    print("-" * 75)

    for xp, yp in points:
        # Find closest grid index
        i = np.argmin(np.abs(y - yp))
        j = np.argmin(np.abs(x - xp))

        numerical = T_num[i, j]
        exact = T_exact[i, j]
        abs_error = abs(numerical - exact)

        print(f"{x[j]:8.3f} {y[i]:8.3f} {numerical:15.8e} {exact:15.8e} {abs_error:15.8e}")


# -------------------------
# Main script
# -------------------------

# User inputs
Lx = 1.0
Ly = 1.0
nx = 41
ny = 41
alpha = 3.352 # thermal diffusivity (steel, stainless 310 at 25 °C)
dt = 1e-5
t_final = 0.01
T0 = 100 # maximum initial temperature at the center of the plate

# Run FTCS solver
x, y, X, Y, T_num, dx, dy, rx, ry = ftcs_2d_heat(
    Lx=Lx,
    Ly=Ly,
    nx=nx,
    ny=ny,
    alpha=alpha,
    dt=dt,
    t_final=t_final,
    T0=T0
)

# Compute analytical solution
T_exact = analytical_solution(X, Y, t_final, alpha, Lx, Ly, T0)

# Compute errors
L2_error, max_error, error = compute_errors(T_num, T_exact)

# Print summary
print("FTCS Explicit Method Results")
print("-" * 40)
print(f"Domain: [0, {Lx}] x [0, {Ly}]")
print(f"Grid size: nx = {nx}, ny = {ny}")
print(f"dx = {dx:.6f}, dy = {dy:.6f}")
print(f"dt = {dt:.6e}")
print(f"Final time = {t_final}")
print(f"rx = {rx:.6f}, ry = {ry:.6f}")
print(f"rx + ry = {rx + ry:.6f}")
print(f"L2 error = {L2_error:.8e}")
print(f"Maximum error = {max_error:.8e}")

# Print comparison table
print_selected_values(x, y, T_num, T_exact)

# Plot results
plot_heatmap(X, Y, T_num, "FTCS Numerical Solution")
plot_heatmap(X, Y, T_exact, "Analytical Solution")
#plot_heatmap(X, Y, np.abs(error), "Absolute Error")
plt.show()