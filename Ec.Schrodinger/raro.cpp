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
    // --- Parámetros de simulación ---
    const int m = 1000;               
    const double k = 0.0014;        // Delta t
    const double a = 0.0;           
    const double b = PI;            
    const double T_final = 1.0;     
    const int N_steps = round(T_final / k); 
    const cd eye(0.0, 1.0);

    double h = (b - a) / static_cast<double>(m); 
    double lam = k / (pow(h, 2));

    cout << "--- Simulación: i*u_t = -u_xx | u = exp(x + it) ---" << endl;
    cout << "Delta x: " << h << " | Delta t: " << k << " | Pasos: " << N_steps << endl;

    auto start_time = chrono::high_resolution_clock::now();

    // 1. Condición Inicial: u(x,0) = e^x
    VectorXcd U(m + 1);
    for (int i = 0; i <= m; ++i) {
        double x = a + i * h;
        U(i) = exp(x); 
    }

    ofstream outfile("raro_1000.dat");
    outfile << "# t            x            u_real       u_imag       u_abs" << endl;
    outfile << fixed << setprecision(6);

    // 2. Matrices A y B (Crank-Nicolson para i*u_t = -u_xx)
    int n_int = m - 1;
    MatrixXcd A = MatrixXcd::Zero(n_int, n_int);
    MatrixXcd B = MatrixXcd::Zero(n_int, n_int);
    
    cd diag_A = 2.0 * eye - lam;
    cd diag_B = 2.0 * eye + lam;

    for (int i = 0; i < n_int; ++i) {
        A(i, i) = diag_A; B(i, i) = diag_B;
        if (i > 0) { 
            A(i, i - 1) = lam / 2.0; 
            B(i, i - 1) = -lam / 2.0; 
        }
        if (i < n_int - 1) { 
            A(i, i + 1) = lam / 2.0; 
            B(i, i + 1) = -lam / 2.0; 
        }
    }
    auto solver = A.partialPivLu();

    // 3. Evolución Temporal
    double t_final_efectivo = 0.0;
    for (int nn = 0; nn < N_steps; ++nn) {
        double t_n = nn * k;
        double t_np1 = (nn + 1) * k;
        t_final_efectivo = t_np1;

        cd bc_L_n = exp(a + eye * t_n);
        cd bc_R_n = exp(b + eye * t_n);
        cd bc_L_np1 = exp(a + eye * t_np1);
        cd bc_R_np1 = exp(b + eye * t_np1);

        VectorXcd C = VectorXcd::Zero(n_int);
        C(0) = -(lam / 2.0) * (bc_L_n + bc_L_np1);
        C(n_int - 1) = -(lam / 2.0) * (bc_R_n + bc_R_np1);

        U.segment(1, n_int) = solver.solve((B * U.segment(1, n_int)) + C);
        U(0) = bc_L_np1;
        U(m) = bc_R_np1;

        for (int i = 0; i <= m; ++i) {
            outfile << setw(13) << t_np1 << setw(13) << a + i * h 
                    << setw(13) << U(i).real() << setw(13) << U(i).imag() 
                    << setw(13) << abs(U(i)) << endl;
        }
        outfile << endl;
    }

    auto end_time = chrono::high_resolution_clock::now();
    
    // --- CÁLCULO DEL ERROR GLOBAL (NORMA L2) ---
    double suma_error_sq = 0.0;
    for (int i = 0; i <= m; ++i) {
        double x_i = a + i * h;
        // Solución analítica para i*u_t = -u_xx con u(x,0)=e^x es u(x,t) = exp(x + i*t)
        cd u_analitica = exp(x_i + eye * t_final_efectivo);
        
        // Sumar el cuadrado de la magnitud de la diferencia
        suma_error_sq += pow(abs(U(i) - u_analitica), 2);
    }
    double error_global = sqrt(suma_error_sq);

    // --- SALIDA DE RESULTADOS ---
    cout << "------------------------------------" << endl;
    cout << "Finalizado en Debian." << endl;
    cout << "Tiempo de ejecución: " << fixed << setprecision(6) 
         << chrono::duration<double>(end_time - start_time).count() << " s" << endl;
    cout << "Error Global (Norma L2): " << scientific << setprecision(6) << error_global << endl;
    cout << "------------------------------------" << endl;
    
    outfile.close();
    return 0;
}