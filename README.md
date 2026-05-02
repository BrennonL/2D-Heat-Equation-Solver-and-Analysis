# 2D Heat Equation Solver and Analysis

This repository contains the implementation and analysis of numerical methods for solving the **two-dimensional transient heat equation**, developed as part of the **AOE 5404: Numerical Methods Course Project**.

---

## 📌 Project Overview

The objective of this project is to study and compare **explicit and implicit finite difference methods** for solving the 2D heat diffusion equation:

$$
\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right)
$$

The project includes:
- Implementation of numerical solvers  
- Verification using an analytical solution  
- Error and stability analysis  
- Visualization of temperature evolution  

---

## ⚙️ Methods Implemented

### Explicit Method
- **FTCS (Forward-Time Central-Space)**
    - Simple and computationally efficient
    - Conditionally stable

### Implicit Methods
- **Backward Euler**
- **Crank–Nicolson**
    - Unconditionally stable (allow larger time steps)
    - Require solving linear systems at each time step

---

## 🔢 Numerical Method Formulations

### 🔹 Explicit Method (FTCS)

The FTCS (Forward-Time Central-Space) scheme is given by:

$$T_{i,j}^{n+1} = T_{i,j}^{n} + r_x \left(T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}\right) + r_y \left(T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}\right)$$


where:

$$
r_x = \frac{\alpha \Delta t}{\Delta x^2}, \qquad r_y = \frac{\alpha \Delta t}{\Delta y^2}
$$

**Stability condition:**

$$
r_x + r_y \le \frac{1}{2}
$$

---

### 🔹 Implicit Method (Backward Euler)

The Backward Euler method uses a fully implicit time discretization:

$$T_{i,j}^{n+1} = T_{i,j}^{n} + r_x \left(T_{i+1,j}^{n+1} - 2T_{i,j}^{n+1} + T_{i-1,j}^{n+1}\right) + r_y \left(T_{i,j+1}^{n+1} - 2T_{i,j}^{n+1} + T_{i,j-1}^{n+1}\right)$$

Matrix form:

$$
A \mathbf{T}^{n+1} = \mathbf{T}^n + \mathbf{b}
$$

**Property:**
- Unconditionally stable

---

### 🔹 Implicit Method (Crank–Nicolson)

The Crank–Nicolson method averages FTCS and Backward Euler:

$$T_{i,j}^{n+1} = T_{i,j}^{n} + \frac{r_x}{2} \left[\left(T_{i+1,j}^{n} - 2T_{i,j}^{n} + T_{i-1,j}^{n}\right) + \left(T_{i+1,j}^{n+1} - 2T_{i,j}^{n+1} + T_{i-1,j}^{n+1}\right)\right] + \frac{r_y}{2} \left[\left(T_{i,j+1}^{n} - 2T_{i,j}^{n} + T_{i,j-1}^{n}\right) + \left(T_{i,j+1}^{n+1} - 2T_{i,j}^{n+1} + T_{i,j-1}^{n+1}\right)\right]$$

Matrix form:

$$
A \mathbf{T}^{n+1} = B \mathbf{T}^{n} + \mathbf{b}
$$

**Properties:**
- Unconditionally stable  
- Second-order accurate in time  
---

## 🧪 Verification

The numerical solutions are validated against the analytical solution:

$$
T(x,y,t) = T_0 e^{-\alpha \pi^2 \left(\frac{1}{L_x^2} + \frac{1}{L_y^2}\right)t}
\sin\left(\frac{\pi x}{L_x}\right)
\sin\left(\frac{\pi y}{L_y}\right)
$$

Error metrics used:
- **L2 norm**
- **Maximum error**

---

## 📊 Features

- 2D heat equation solver (explicit & implicit)
- Stability and accuracy analysis
- Heatmap visualization of temperature fields
- Modular Python implementation using NumPy and Matplotlib

---
## 👥 Contributors
Brennon Laney  
Hatice Nur Sun  
Khoa Nguyen  

---
## 📚 Course
AOE 5404 — Numerical Methods in Aerospace Engineering