import numpy as np
import matplotlib.pyplot as plt

def exact_solution(X, Y, t, alpha, Lx, Ly):
    return (
        np.sin(np.pi * X / Lx)
        * np.sin(np.pi * Y / Ly)
        * np.exp(-alpha * np.pi**2 * (1 / Lx**2 + 1 / Ly**2) * t)
    )

def main():
  ...

if __name__ == "__main__":
  main()