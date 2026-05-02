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

def exact_solution(X, Y, t, alpha, Lx, Ly):
    return (
        np.sin(np.pi * X / Lx)
        * np.sin(np.pi * Y / Ly)
        * np.exp(-alpha * np.pi**2 * (1 / Lx**2 + 1 / Ly**2) * t)
    )

def main():
  Lx = 1.0
  Ly = 1.0
  nx = 41
  ny = 41
  alpha = 3.352 # thermal diffusivity (steel, stainless 310 at 25 °C)
  dt = 1e-5
  t_final = 0.01
  T0 = 100 # maximum initial temperature at the center of the plate
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
  b_euler = backward_Euler(dx,dt,t_final,ALPHA,T_flat, (0,Lx,0,Ly),(0,0,0,0))
  plt.clf()
  plt.imshow(b_euler, origin='lower', cmap='hot')
  plt.colorbar(label="Temperature")
  plt.title("Heatmap")
  plt.xlabel("x")
  plt.ylabel("y")
  plt.savefig("backward-euler.png")
  plt.clf()



if __name__ == "__main__":
  main()