import numpy as np
import math
import pandas as pd
from scipy.linalg import solve_banded
import time
import os

def solve_heat_equation_cn(L, T_final, alpha, N, M, f, bc_left=0.0, bc_right=0.0):
    """
    Solves the 1D Heat Equation using the Crank-Nicolson method.

    Args:
        L: Length of the rod.
        T_final: Final time.
        alpha: Thermal diffusivity.
        N: Number of spatial divisions.
        M: Number of time steps.
        f: Initial condition function f(x).
        bc_left: Dirichlet boundary condition at x=0.
        bc_right: Dirichlet boundary condition at x=L.

    Returns:
        x: Spatial grid.
        t: Time grid.
        U: Solution matrix (N+1, M+1).
    """
    dx = L / N
    dt = T_final / M
    r = alpha**2 * dt / (dx**2)

    x = np.linspace(0, L, N+1)
    t = np.linspace(0, T_final, M+1)

    U = np.zeros((N+1, M+1))

    # Initial Condition
    U[:, 0] = [f(xi) for xi in x]

    # Boundary Conditions (Dirichlet)
    U[0, :] = bc_left
    U[N, :] = bc_right

    # Crank-Nicolson matrices (tridiagonal)
    # solve_banded expects (l, u) where l is number of lower diagonals, u is upper
    # Matrix A (Implicit part): (1 + r) on diagonal, -r/2 on off-diagonals
    # Matrix B (Explicit part): (1 - r) on diagonal, r/2 on off-diagonals

    # For N-1 internal nodes
    diag_A = (1 + r) * np.ones(N-1)
    upper_A = -r/2 * np.ones(N-1)
    lower_A = -r/2 * np.ones(N-1)

    # Banded matrix format for solve_banded: (upper, middle, lower)
    ab = np.zeros((3, N-1))
    ab[0, 1:] = upper_A[:-1]
    ab[1, :] = diag_A
    ab[2, :-1] = lower_A[1:]

    for j in range(M):
        # Calculate RHS = B * U[1:N, j] + boundary contributions
        rhs = np.zeros(N-1)
        for i in range(1, N):
            # i-1 because rhs is for internal nodes 1 to N-1
            idx = i - 1
            # Current internal nodes
            u_prev = U[i-1, j]
            u_curr = U[i, j]
            u_next = U[i+1, j]

            rhs[idx] = (r/2) * u_prev + (1 - r) * u_curr + (r/2) * u_next

        # Add boundary terms for the next time step (j+1)
        rhs[0] += (r/2) * U[0, j+1]
        rhs[-1] += (r/2) * U[N, j+1]

        # Solve A * U[1:N, j+1] = rhs
        U[1:N, j+1] = solve_banded((1, 1), ab, rhs)

    return x, t, U

if __name__ == "__main__":
    # Parameters
    L = 1.0
    T_final = 0.1
    alpha = 1.0
    N = 50
    M = 100

    def initial_condition(x):
        return math.sin(math.pi * x)

    print(f"Solving 1D Heat Equation...")
    print(f"L={L}, T_final={T_final}, alpha={alpha}")
    print(f"N={N} (spatial), M={M} (temporal)")

    start_time = time.time()
    x_grid, t_grid, U_sol = solve_heat_equation_cn(L, T_final, alpha, N, M, initial_condition)
    end_time = time.time()

    print(f"Execution time: {end_time - start_time:.6f} seconds")

    # Analytical solution for comparison: u(x,t) = exp(-(alpha*pi)^2 * t) * sin(pi * x)
    def analytical(x, t, alpha_val):
        return math.exp(-((alpha_val * math.pi)**2) * t) * math.sin(math.pi * x)

    # Compare at T_final
    errors = []
    for i, xi in enumerate(x_grid):
        approx = U_sol[i, -1]
        exact = analytical(xi, T_final, alpha)
        errors.append(abs(approx - exact))

    max_error = max(errors)
    print(f"Max error at T_final: {max_error:.2e}")

    # Save results
    output_filename = "results_python.csv"
    df = pd.DataFrame(U_sol.T, columns=[f"x_{xi:.2f}" for xi in x_grid])
    df.insert(0, "time", t_grid)
    df.to_csv(output_filename, index=False)
    print(f"Results saved to {output_filename}")
