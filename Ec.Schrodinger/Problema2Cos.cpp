#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <chrono> // Librería para medir tiempos

using namespace std;
using namespace Eigen;

typedef complex<double> cd;
typedef Matrix<cd, Dynamic, 1> VectorXcd;
typedef Matrix<cd, Dynamic, Dynamic> MatrixXcd;

// Definición de pi para mayor precisión
const double PI = acos(-1.0);

int main() {
    // --- Parámetros del Problema 2 ---
    const int m = 1000;              
    const double k = 0.0011;        // Delta t
    const double a = -PI;           // Límite izquierdo
    const double b = PI;            // Límite derecho
    const double T_final = 1.0;     // Tiempo final simulado
    const int N_steps = round(T_final / k); 
    const cd eye(0.0, 1.0);

    // Delta x calculado a partir de m
    double h = (b - a) / static_cast<double>(m); 
    double lam = k / (pow(h, 2));

    cout << "--- Configuración Problema 2 ---" << endl;
    cout << "Delta x (h): " << h << " | Delta t (k): " << k << endl;
    cout << "Lambda: " << lam << " | Pasos totales: " << N_steps << endl;

    // Medición de tiempo de inicio
    auto start_time = chrono::high_resolution_clock::now();

    // 1. Inicialización de la función de onda (U completo)
    // Condición inicial: u(x,0) = cos(x)
    VectorXcd U(m + 1);
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        U(i) = cos(x); 
    }

    // Preparar archivo de salida
    ofstream outfile("resultados_p2_m1000.dat");
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

    // 2. Construcción de matrices A y B (Sistema tridiagonal)
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

    // Pre-factorización para optimizar el bucle
    auto solver = A.partialPivLu();

    // 3. Evolución Temporal
    double t_actual_val = 0.0; 
    for (int nn = 0; nn < N_steps; ++nn) {
        double t_n = nn * k;
        double t_np1 = (nn + 1) * k;
        t_actual_val = t_np1;

        // Fronteras basadas en u(x,t) = exp(it) * cos(x)
        cd bc_L_n = exp(eye * t_n) * cos(a);
        cd bc_R_n = exp(eye * t_n) * cos(b);
        cd bc_L_np1 = exp(eye * t_np1) * cos(a);
        cd bc_R_np1 = exp(eye * t_np1) * cos(b);

        VectorXcd C = VectorXcd::Zero(n_int);
        C(0) = lam * (bc_L_n + bc_L_np1);
        C(n_int - 1) = lam * (bc_R_n + bc_R_np1);

        // Resolver el paso de tiempo
        VectorXcd rhs = (B * U.segment(1, n_int)) + C;
        VectorXcd sol = solver.solve(rhs);

        // Actualizar vector de solución
        U(0) = bc_L_np1;
        U(m) = bc_R_np1;
        U.segment(1, n_int) = sol;

        // Guardar en el .dat
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

    // --- CÁLCULO DEL ERROR GLOBAL (NORMA L2) ---
    double suma_error_sq = 0.0;
    for (int i = 0; i <= m; ++i) {
        double x_i = a + i * h;
        // Solución analítica: u(x, t) = exp(i*t) * cos(x)
        cd u_analitica = exp(eye * t_actual_val) * cos(x_i);
        
        // Diferencia al cuadrado (magnitud compleja)
        suma_error_sq += pow(abs(U(i) - u_analitica), 2);
    }
    double error_global = sqrt(suma_error_sq);

    // Medición de tiempo final
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    outfile.close();

    cout << "\n--- Finalizado ---" << endl;
    cout << "Archivo 'resultados_p2_m1000.dat' generado." << endl;
    cout << "Tiempo total: " << fixed << setprecision(6) << elapsed.count() << " segundos." << endl;
    cout << "Error Global (Norma L2): " << scientific << setprecision(6) << error_global << endl;

    return 0;
}