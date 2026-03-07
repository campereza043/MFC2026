#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <chrono> 

using namespace std;
using namespace Eigen;

typedef complex<double> cd;
typedef Matrix<cd, Dynamic, 1> VectorXcd;
typedef Matrix<cd, Dynamic, Dynamic> MatrixXcd;

const double PI = acos(-1.0);

int main() {
    // --- Parámetros ajustados ---
    const int m = 1000;               
    const double k = 0.0014;        // Delta t
    const double a = 0.0;           
    const double b = PI;            
    const double T_final = 1.0;     
    const int N_steps = round(T_final / k); 
    const cd eye(0.0, 1.0);

    double h = (b - a) / static_cast<double>(m); 
    double lam = k / (pow(h, 2));

    cout << "--- Configuración Problema 3 (u = e^{it+x}) ---" << endl;
    cout << "Delta x (h): " << h << " | Delta t (k): " << k << endl;
    cout << "Pasos totales: " << N_steps << endl;

    auto start_time = chrono::high_resolution_clock::now();

    // 1. Condición Inicial: u(x,0) = e^{x}
    VectorXcd U(m + 1);
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        U(i) = exp(x); 
    }

    ofstream outfile("resultados_p3_m1000.dat");
    outfile << "# t            x            u_real       u_imag       u_abs" << endl;
    outfile << fixed << setprecision(6);

    // Guardar t = 0
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        outfile << setw(13) << 0.0 << setw(13) << x << setw(13) << U(i).real() 
                << setw(13) << U(i).imag() << setw(13) << abs(U(i)) << endl;
    }
    outfile << endl;

    // 2. Matrices A y B
    int n_int = m - 1;
    MatrixXcd A = MatrixXcd::Zero(n_int, n_int);
    MatrixXcd B = MatrixXcd::Zero(n_int, n_int);
    cd diag_A = 2.0 * eye + 2.0 * lam;
    cd diag_B = 2.0 * eye - 2.0 * lam;

    for (int i = 0; i < n_int; ++i) {
        A(i, i) = diag_A; B(i, i) = diag_B;
        if (i > 0) { A(i, i - 1) = -lam; B(i, i - 1) = lam; }
        if (i < n_int - 1) { A(i, i + 1) = -lam; B(i, i + 1) = lam; }
    }
    auto solver = A.partialPivLu();

    // 3. Evolución Temporal
    for (int nn = 0; nn < N_steps; ++nn) {
        double t_n = nn * k;
        double t_np1 = (nn + 1) * k;

        // Fronteras basadas en u(x,t) = exp(i*t + x)
        cd bc_L_n = exp(eye * t_n + a);
        cd bc_R_n = exp(eye * t_n + b);
        cd bc_L_np1 = exp(eye * t_np1 + a);
        cd bc_R_np1 = exp(eye * t_np1 + b);

        VectorXcd C = VectorXcd::Zero(n_int);
        C(0) = lam * (bc_L_n + bc_L_np1);
        C(n_int - 1) = lam * (bc_R_n + bc_R_np1);

        U.segment(1, n_int) = solver.solve((B * U.segment(1, n_int)) + C);
        U(0) = bc_L_np1;
        U(m) = bc_R_np1;

        for (int i = 0; i <= m; ++i) {
            outfile << setw(13) << t_np1 << setw(13) << a + i * h << setw(13) << U(i).real() 
                    << setw(13) << U(i).imag() << setw(13) << abs(U(i)) << endl;
        }
        outfile << endl;
    }

    auto end_time = chrono::high_resolution_clock::now();
    cout << "Tiempo de ejecución: " << chrono::duration<double>(end_time - start_time).count() << " s" << endl;
    outfile.close();
    return 0;
}