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
    // --- Parámetros del Problema 3 ---
    const int m = 20;               // M = 20
    const double k = 0.0014;        // Delta t
    const double a = 0.0;           // Límite izquierdo
    const double b = PI;            // Límite derecho
    const double T_final = 1.0;     // Tiempo simulado (puedes ajustarlo)
    const int N_steps = round(T_final / k); 
    const cd eye(0.0, 1.0);

    // Delta x (h) será aproximadamente 0.1571 para M=20 en [0, PI]
    double h = (b - a) / static_cast<double>(m); 
    double lam = k / (pow(h, 2));

    cout << "--- Configuración Problema 3 ---" << endl;
    cout << "Delta x (h): " << h << " | Delta t (k): " << k << endl;
    cout << "Lambda: " << lam << " | Pasos totales: " << N_steps << endl;

    auto start_time = chrono::high_resolution_clock::now();

    // 1. Inicialización de la función de onda
    // Condición inicial del Problema 3: u(x,0) = x
    VectorXcd U(m + 1);
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        U(i) = cd(x, 0.0); 
    }

    // Preparar archivo de salida
    ofstream outfile("eitx.dat");
    outfile << "# t            x            u_real       u_imag       u_abs" << endl;
    outfile << fixed << setprecision(6);

    // Guardar t = 0
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        outfile << setw(13) << 0.0 
                << setw(13) << x 
                << setw(13) << U(i).real() 
                << setw(13) << U(i).imag() 
                << setw(13) << abs(U(i)) << endl;
    }
    outfile << endl;

    // 2. Construcción de matrices A y B
    int n_int = m - 1;
    MatrixXcd A = MatrixXcd::Zero(n_int, n_int);
    MatrixXcd B = MatrixXcd::Zero(n_int, n_int);

    cd diag_A = 2.0 * eye + 2.0 * lam;
    cd diag_B = 2.0 * eye - 2.0 * lam;

    for (int i = 0; i < n_int; ++i) {
        A(i, i) = diag_A;
        B(i, i) = diag_B;
        if (i > 0) {
            A(i, i - 1) = -lam;
            B(i, i - 1) = lam;
        }
        if (i < n_int - 1) {
            A(i, i + 1) = -lam;
            B(i, i + 1) = lam;
        }
    }

    auto solver = A.partialPivLu();

    // 3. Evolución Temporal
    for (int nn = 0; nn < N_steps; ++nn) {
        double t_actual = nn * k;
        double t_np1 = (nn + 1) * k;

        // Fronteras basadas en u(x,t) = exp(it) * x
        cd bc_L_n = exp(eye * t_actual) * a;
        cd bc_R_n = exp(eye * t_actual) * b;
        cd bc_L_np1 = exp(eye * t_np1) * a;
        cd bc_R_np1 = exp(eye * t_np1) * b;

        VectorXcd C = VectorXcd::Zero(n_int);
        C(0) = lam * (bc_L_n + bc_L_np1);
        C(n_int - 1) = lam * (bc_R_n + bc_R_np1);

        VectorXcd rhs = (B * U.segment(1, n_int)) + C;
        VectorXcd sol = solver.solve(rhs);

        U(0) = bc_L_np1;
        U(m) = bc_R_np1;
        U.segment(1, n_int) = sol;

        // Guardar cada paso de tiempo
        for (int i = 0; i <= m; ++i) {
            double x = a + i * h;
            outfile << setw(13) << t_np1 
                    << setw(13) << x 
                    << setw(13) << U(i).real() 
                    << setw(13) << U(i).imag() 
                    << setw(13) << abs(U(i)) << endl;
        }
        outfile << endl;
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    outfile.close();

    cout << "\n--- Finalizado ---" << endl;
    cout << "Archivo 'eitx.dat' generado exitosamente." << endl;
    cout << "Tiempo de cálculo y E/S: " << elapsed.count() << " segundos." << endl;

    return 0;
}