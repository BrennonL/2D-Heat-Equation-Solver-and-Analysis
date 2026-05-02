from typing import Callable
from numpy import (
  array,ndarray,sin,cos,zeros,fill_diagonal
)
from numpy.linalg import solve
from scipy.optimize import (
  fsolve
)
# from sympy import (
#   symbols,Eq,solve,I
# )

def backward_Euler(
    h:float|int, # detla x and delta y
    delta_t:float|int,
    t:float|int,
    c:float,
    heat_map:array[float,float,float,float],
    boundary_conditions:array[float,float,float,float]
):
  """
  The Backward Euler algorithm for 2D heatmap
  :param h: The stepsize of x and y
  :param delta_t: t stepsize
  :param t: The simulation time duration
  :param c: The C constant (specific to different materials)
  :param heat_map: The dimensions of the 2D plate
  :param boundary_conditions: The heat of the boundary
  :returns: The final heat map after t time
  """
  (x_min, x_max, y_min, y_max) = heat_map # The size of the plate
  (Tl, Tr, Tt, Tb) = boundary_conditions # The temperature of the edges
  Lambda = (c*delta_t) / (h**2) #! This only works if delta x and delta y are the same
  nx = int((x_max - x_min) / h -1)
  ny = int((y_max - y_min) / h -1)
  N = nx * ny

  initial_x:ndarray = zeros(int(N+1))
  # TODO -- Implement boundary conditions
  initial_x[0] = ...

  # Create sparce matrix or coefficient matrix
  A:ndarray = zeros((int(N+1),int(N+1)))
  b:ndarray = zeros(int(N+1))
  fill_diagonal(A, (1 + 2*Lambda + 2*Lambda))
  for i in range(len(A)): # rows
    for j in range(len(A[i])): # columns
      if A[i][j] == 1 + 2*Lambda:
        if j>0: # Left
          A[i][j-1] = - Lambda
          b[i] += Tl
        if j<(len(A[i]) - 1): # right
          A[i][j+1] = - Lambda
          b[i] += Tr
        if j-nx>0: # bottom
          A[i][j-nx] = - Lambda
          b[i] += Tb
        if j+nx<(len(A[i])-1): # top
          A[i][j+nx] = - Lambda
          b[i] += Tt

  # Run simulation
  n = t / delta_t # Number of time steps
  curr_x = initial_x
  for k in range(int(n)):
    next_x = solve(A, b)
    next_x[0] = 100
    curr_x = next_x
  return curr_x

def crank_nicolson(
    h:float, # detla x and delta y
    delta_t,
    t,
    c:float,
    heat_map:array[float,float,float,float],
    boundary_conditions:array[float,float,float,float]
)->ndarray:
  """
  The Crank Nicolson algorithm for 2D heatmap
  :param h: The stepsize of x and y
  :param delta_t: t stepsize
  :param t: The simulation time duration
  :param c: The C constant (specific to different materials)
  :param heat_map: The dimensions of the 2D plate
  :param boundary_conditions: The heat of the boundary
  :returns: The final heat map after t time
  """
  (x_min, x_max, y_min, y_max) = heat_map # The size of the plate
  (Tl, Tr, Tt, Tb) = boundary_conditions # The temperature of the edges
  Lambda = (c*delta_t) / (h**2)
  nx = int((x_max - x_min) / h -1)
  ny = int((y_max - y_min) / h -1)
  N = nx * ny

  initial_x:ndarray = zeros(int(N+1))
  # TODO -- Implement boundary conditions
  initial_x[0] = ...

  # Create sparce matrix or coefficient matrix
  A:ndarray = zeros((int(N+1),int(N+1)))
  b:ndarray = zeros(int(N+1))
  fill_diagonal(A, (1 + 2*Lambda))
  for i in range(len(A)): # rows
    for j in range(len(A[i])): # columns
      if A[i][j] == 1 + 2*Lambda:
        if j>0: # Left
          A[i][j-1] = - Lambda
          b[i] += Tl
        if j<(len(A[i]) - 1): # right
          A[i][j+1] = - Lambda
          b[i] += Tr
        if j-nx>0: # bottom
          A[i][j-nx] = - Lambda
          b[i] += Tb
        if j+nx<(len(A[i])-1): # top
          A[i][j+nx] = - Lambda
          b[i] += Tt

  # Run simulation
  n = t / delta_t # Number of time steps
  curr_x = initial_x
  for k in range(int(n)):
    next_x = solve(A, b)
    next_x[0] = 100
    curr_x = next_x
  return curr_x