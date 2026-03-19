#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <chrono>

using namespace std;

/**
 * Solves a tridiagonal system Ax = d using the Thomas Algorithm.
 *
 * @param a Subdiagonal elements (a[0] not used).
 * @param b Main diagonal elements.
 * @param c Superdiagonal elements (c[n-1] not used).
 * @param d RHS vector.
 * @return Solution vector x.
 */
vector<double> thomas_algorithm(const vector<double>& a, const vector<double>& b, const vector<double>& c, vector<double> d) {
    int n = d.size();
    vector<double> c_prime(n, 0.0);
    vector<double> d_prime(n, 0.0);
    vector<double> x(n, 0.0);

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double m = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / m;
    }

    // Back substitution
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}

int main() {
    // Parameters
    const double L = 1.0;
    const double T_final = 0.1;
    const double alpha = 1.0;
    const int N = 50;       // Spatial steps
    const int M = 100;      // Time steps
    const double dx = L / N;
    const double dt = T_final / M;
    const double r = (alpha * alpha * dt) / (dx * dx);

    cout << "Solving 1D Heat Equation (Crank-Nicolson) in C++" << endl;
    cout << "L=" << L << ", T_final=" << T_final << ", alpha=" << alpha << endl;
    cout << "N=" << N << ", M=" << M << ", r=" << r << endl;

    // Initial and boundary conditions
    vector<double> U(N + 1);
    for (int i = 0; i <= N; i++) {
        double x = i * dx;
        U[i] = sin(M_PI * x); // f(x) = sin(pi*x)
    }

    // Boundary conditions: Dirichlet (fixed at 0)
    U[0] = 0.0;
    U[N] = 0.0;

    // Tridiagonal system for internal nodes (1 to N-1)
    int n_int = N - 1;
    vector<double> a(n_int, -r / 2.0);
    vector<double> b(n_int, 1.0 + r);
    vector<double> c(n_int, -r / 2.0);
    vector<double> d(n_int);

    string output_filename = "results_cpp.dat";
    ofstream outfile(output_filename);
    outfile << fixed << setprecision(6);

    auto start = chrono::high_resolution_clock::now();

    for (int j = 0; j < M; j++) {
        // Build RHS for internal nodes
        for (int i = 1; i < N; i++) {
            int idx = i - 1;
            d[idx] = (r / 2.0) * U[i - 1] + (1.0 - r) * U[i] + (r / 2.0) * U[i + 1];
        }

        // Solve for next time step
        vector<double> U_next_int = thomas_algorithm(a, b, c, d);

        // Update U (internal nodes)
        for (int i = 1; i < N; i++) {
            U[i] = U_next_int[i - 1];
        }

        // Save current time step to file
        double t_curr = (j + 1) * dt;
        outfile << t_curr << " ";
        for (int i = 0; i <= N; i++) {
            outfile << U[i] << (i == N ? "" : " ");
        }
        outfile << endl;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;

    cout << "Execution time: " << diff.count() << " seconds" << endl;
    outfile.close();

    // Verification with analytical solution at T_final: exp(-(alpha*pi)^2 * T_final) * sin(pi * x)
    double max_err = 0.0;
    for (int i = 0; i <= N; i++) {
        double x = i * dx;
        double exact = exp(-pow(alpha * M_PI, 2) * T_final) * sin(M_PI * x);
        double err = abs(U[i] - exact);
        if (err > max_err) max_err = err;
    }
    cout << "Max error at T_final: " << max_err << endl;
    cout << "Results saved to " << output_filename << endl;

    return 0;
}
