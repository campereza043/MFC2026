#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

typedef complex<double> cd;
typedef Matrix<cd, Dynamic, 1> VectorXcd;
typedef Matrix<cd, Dynamic, Dynamic> MatrixXcd;

int main() {
    // Parámetros idénticos a tu código de Fortran
    const int m = 20;
    const double k = 0.0013;
    const double a = -2.0;
    const double b = 2.0;
    const int N_steps = 769; 
    const cd eye(0.0, 1.0);

    double h = (b - a) / static_cast<double>(m);
    double lam = k / (pow(h, 2));

    // Inicialización de la función de onda
    VectorXcd U(m + 1);
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        U(i) = sin(x);
    }

    // Preparar archivo de salida
    ofstream outfile("resultados_cpp.dat");
    outfile << "# t            x            u_real       u_imag       u_abs" << endl;
    outfile << fixed << setprecision(6);

    // Guardar condición inicial (t = 0)
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        outfile << setw(13) << 0.0 
                << setw(13) << x 
                << setw(13) << U(i).real() 
                << setw(13) << U(i).imag() 
                << setw(13) << abs(U(i)) << endl;
    }
    outfile << endl;

    // Construcción de matrices A y B (Tridiagonales)
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

    // Evolución Temporal (Guardando todos los pasos como en Fortran)
    for (int nn = 0; nn < N_steps; ++nn) {
        double t_actual = nn * k;
        double t_np1 = (nn + 1) * k;

        // Condiciones de frontera
        cd bc_L_n = exp(eye * t_actual) * sin(a);
        cd bc_R_n = exp(eye * t_actual) * sin(b);
        cd bc_L_np1 = exp(eye * t_np1) * sin(a);
        cd bc_R_np1 = exp(eye * t_np1) * sin(b);

        // Vector de carga C
        VectorXcd C = VectorXcd::Zero(n_int);
        C(0) = lam * (bc_L_n + bc_L_np1);
        C(n_int - 1) = lam * (bc_R_n + bc_R_np1);

        // Resolver sistema
        VectorXcd U_inner = U.segment(1, n_int);
        VectorXcd rhs = (B * U_inner) + C;
        VectorXcd sol = solver.solve(rhs);

        // Actualizar U
        U(0) = bc_L_np1;
        U(m) = bc_R_np1;
        U.segment(1, n_int) = sol;

        // Guardar cada paso de tiempo en el archivo .dat
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

    outfile.close();
    cout << "¡Archivo resultados_cpp.dat generado con los mismos pasos que Fortran!" << endl;

    return 0;
}