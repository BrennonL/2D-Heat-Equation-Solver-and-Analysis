import numpy as np
import matplotlib.pyplot as plt
from explicit_methods import (
  analytical_solution,
  ftcs_2d_heat
)
from implicit_methods import (
   crank_nicolson,
   backward_Euler
)
from constants import ALPHA

def main():
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
      alpha=ALPHA,
      dt=dt,
      t_final=t_final,
      T0=T0
  )
  
  # Compute analytical solution
  T_exact = analytical_solution(X, Y, t_final, alpha, Lx, Ly, T0)

  # Grid setup
  dx = Lx / (nx - 1) 
  dy = Ly / (ny - 1)
  x = np.linspace(0, Lx, nx)
  y = np.linspace(0, Ly, ny)
  X, Y = np.meshgrid(x, y)

  # Stability parameters
  rx = alpha * dt / dx**2
  ry = alpha * dt / dy**2

  # Initial condition
  T = T0 * np.sin(np.pi * X / Lx) * np.sin(np.pi * Y / Ly)
  T = T[1:-1, 1:-1]
  T_flat = T.flatten(order="F")

  # Run Backward Euler
  b_euler = backward_Euler(
      dx,
      1e-5,
      t_final,
      alpha,
      T_flat,
      (0, Lx, 0, Ly),
      (0, 0, 0, 0)
  )

  # Put interior solution back into full grid with zero boundaries
  T_be_full = np.zeros((ny, nx))
  T_be_full[1:-1, 1:-1] = b_euler

  # Compare max temperature between methods
  print("FTCS max:", np.max(T_num))
  print("Backward Euler max:", np.max(T_be_full))
  print("Exact max:", np.max(T_exact))

  # Plot Backward Euler using contour map
  vmin = 0
  vmax = 55
  plt.figure(figsize=(6, 5))
  plt.contourf(X, Y, T_be_full, levels=30, cmap="hot", vmin=vmin, vmax=vmax)
  plt.colorbar(label="Temperature")
  plt.xlabel("x")
  plt.ylabel("y")
  plt.title("Backward Euler Numerical Solution")
  plt.tight_layout()
  plt.savefig("Backward-Euler-Solution.png")
  plt.show()


if __name__ == "__main__":
  main()