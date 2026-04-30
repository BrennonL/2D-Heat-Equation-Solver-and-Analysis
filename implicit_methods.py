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

def FTCS():
  ...

def backward_Euler(
    x_curr:ndarray,h:float,t:float,f:Callable[[array],array]
):
  """
  Backward Euler Method
  :param x_curr: The current state
  :param h: The step size
  :param t: The current time
  :param f: The fucntion
  :returns: The next state
  """
  coeficient_maxtri:ndarray = array([
    
  ])
  t_next = t + h
  def equation(x_next:ndarray):
    return - x_next + x_curr + h*f(t_next,x_next)
  x_next = fsolve(equation, x_curr)
  return x_next

def crank_nicolson(
    delta_t,
    t,
    delta_x,
    x_max,
    c
):
  Lambda = (c*delta_t) / (delta_x**2)
  N = x_max/delta_x # Matrix size

  initial_x:ndarray = zeros(int(N+1))
  # TODO -- Implement boundary conditions
  initial_x[0] = ...

  # Create sparce matrix or coefficient matrix
  A:ndarray = zeros((int(N+1),int(N+1)))
  fill_diagonal(A, (1 + 2*Lambda))
  for i in range(len(A)): # rows
    for j in range(len(A[i])): # columns
      if A[i][j] == 1 + 2*Lambda:
        if j>0:
          A[i][j-1] = - Lambda
        if j<(len(A[i]) - 1):
          A[i][j+1] = - Lambda

  # Run simulation
  n = t / delta_t # Number of time steps
  curr_x = initial_x
  for k in range(int(n)):
    next_x = solve(A, curr_x)
    next_x[0] = 100
    curr_x = next_x
  return curr_x