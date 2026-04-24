from typing import Callable
from numpy import (
  array,ndarray,sin,cos
)
from scipy.optimize import (
  fsolve
)
from sympy import (
  symbols,Eq,solve,I
)

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
  t_next = t + h
  def equation(x_next:ndarray):
    return - x_next + x_curr + h*f(t_next,x_next)
  x_next = fsolve(equation, x_curr)
  return x_next

def crank_nicolson():
  ...