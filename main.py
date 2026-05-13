import numpy as np
import time
import matplotlib.pyplot as plt
from explicit_methods import (
  analytical_solution,
  compute_errors,
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
  start_time = time.perf_counter()
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
  ftcs_time = time.perf_counter() - start_time
  
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
  start_time = time.perf_counter()
  b_euler = backward_Euler(
      dx,
      1e-5,
      t_final,
      alpha,
      T_flat,
      (0, Lx, 0, Ly),
      (0, 0, 0, 0)
  )
  backward_euler_time = time.perf_counter() - start_time

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
  # plt.show()

  # Run Crank Nicolson
  start_time = time.perf_counter()
  cn = crank_nicolson(
      dx,
      1e-5,
      t_final,
      alpha,
      T_flat,
      (0, Lx, 0, Ly),
      (0, 0, 0, 0)
  )
  crank_nicolson_time = time.perf_counter() - start_time

  # Put interior solution back into full grid with zero boundaries
  T_cn_full = np.zeros((ny, nx))
  T_cn_full[1:-1, 1:-1] = cn

  # Compare max temperature between methods
  print("FTCS max:", np.max(T_num))
  print("Crank Nicolson max:", np.max(T_cn_full))
  print("Exact max:", np.max(T_exact))

  # Plot Crank Nicolson using contour map
  vmin = 0
  vmax = 55
  plt.figure(figsize=(6, 5))
  plt.contourf(X, Y, T_cn_full, levels=30, cmap="hot", vmin=vmin, vmax=vmax)
  plt.colorbar(label="Temperature")
  plt.xlabel("x")
  plt.ylabel("y")
  plt.title("Crank Nicolson Numerical Solution")
  plt.tight_layout()
  plt.savefig("Crank-Nicolson-Solution.png")
  # plt.show()
  print("Execution Time Summary")
  print("-" * 40)
  print(f"FTCS execution time: {ftcs_time:.6e} seconds")
  print(f"Backward Euler execution time: {backward_euler_time:.6e} seconds")
  print(f"Crank Nicolson execution time: {crank_nicolson_time:.6e} seconds")
  print("-" * 40)

 # Convergence study: L2 error vs h
  grid_sizes = [21, 41, 81, 161]

  h_values = []
  l2_errors = []

  for n in grid_sizes:
      nx_conv = n
      ny_conv = n

      dx_conv = Lx / (nx_conv - 1)

      # Choose stable dt and make sure t_final is reached exactly
      dt_guess = 0.1 * dx_conv**2 / alpha
      n_steps_conv = int(np.ceil(t_final / dt_guess))
      dt_conv = t_final / n_steps_conv

      x_conv, y_conv, X_conv, Y_conv, T_conv, dx_conv, dy_conv, rx_conv, ry_conv = ftcs_2d_heat(
          Lx=Lx,
          Ly=Ly,
          nx=nx_conv,
          ny=ny_conv,
          alpha=alpha,
          dt=dt_conv,
          t_final=t_final,
          T0=T0
      )

      T_exact_conv = analytical_solution(
          X_conv,
          Y_conv,
          t_final,
          alpha,
          Lx,
          Ly,
          T0
      )

      L2_conv, max_conv, error_conv = compute_errors(T_conv, T_exact_conv)

      h_values.append(dx_conv)
      l2_errors.append(L2_conv)

      print(f"Grid size = {nx_conv} x {ny_conv}")
      print(f"h = {dx_conv:.6e}")
      print(f"dt = {dt_conv:.6e}")
      print(f"L2 error = {L2_conv:.8e}")
      print("-" * 40)

  h_values = np.array(h_values)
  l2_errors = np.array(l2_errors)

  sort_idx = np.argsort(h_values)
  h_values = h_values[sort_idx]
  l2_errors = l2_errors[sort_idx]

  plt.figure(figsize=(6, 5))
  plt.loglog(h_values, l2_errors, marker="o")
  plt.xlabel("h")
  plt.ylabel("L2 Error")
  plt.title("Convergence Study: L2 Error vs h")
  plt.grid(True, which="both")
  plt.tight_layout()
  plt.savefig("L2-Convergence-Study.png")
  
if __name__ == "__main__":
  main()