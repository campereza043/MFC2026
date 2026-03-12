#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include <chrono>
#include <fstream> 
#include <iomanip> 

using namespace std;
using namespace Eigen;

typedef complex<double> cd;
typedef Matrix<cd, Dynamic, 1> VectorXcd;
typedef Matrix<cd, Dynamic, Dynamic> MatrixXcd;

int main() {
    // Parámetros
    double a = -2.0;
    double b = 2.0;
    double T = 1.0;
    double k = 0.0013; 
    int M = 1000; 
    
    double h = (b - a) / (M + 1);
    int N = round(T / k);
    double lam = k / (pow(h, 2));

    auto start = chrono::high_resolution_clock::now();

    // Condición inicial
    VectorXcd U(M + 2);
    for (int i = 0; i < M + 2; ++i) {
        U(i) = sin(a + i * h);
    }

    // Matrices tridiagonales
    MatrixXcd A = MatrixXcd::Zero(M, M);
    MatrixXcd B = MatrixXcd::Zero(M, M);
    cd diag_A(2.0 * lam, 2.0); 
    cd diag_B(-2.0 * lam, 2.0); 
    cd off_A(-lam, 0.0);
    cd off_B(lam, 0.0);

    for (int i = 0; i < M; ++i) {
        A(i, i) = diag_A;
        B(i, i) = diag_B;
        if (i > 0) { A(i, i - 1) = off_A; B(i, i - 1) = off_B; }
        if (i < M - 1) { A(i, i + 1) = off_A; B(i, i + 1) = off_B; }
    }

    auto solver = A.partialPivLu();

    // Bucle temporal
    for (int n = 0; n < N; ++n) {
        double t_n = n * k;
        double t_np1 = (n + 1) * k;

        cd gL_n = exp(cd(0, 1) * t_n) * sin(a);
        cd gR_n = exp(cd(0, 1) * t_n) * sin(b);
        cd gL_np1 = exp(cd(0, 1) * t_np1) * sin(a);
        cd gR_np1 = exp(cd(0, 1) * t_np1) * sin(b);

        VectorXcd C = VectorXcd::Zero(M);
        C(0) = lam * (gL_n + gL_np1);
        C(M - 1) = lam * (gR_n + gR_np1);

        VectorXcd rhs = (B * U.segment(1, M)) + C;
        U.segment(1, M) = solver.solve(rhs);
        U(0) = gL_np1;
        U(M + 1) = gR_np1;
    }

    auto end = chrono::high_resolution_clock::now();
    
    // --- CÁLCULO DEL ERROR GLOBAL (L2) ---
    double suma_error_sq = 0.0;
    for (int i = 0; i < M + 2; ++i) {
        double x_i = a + i * h;
        // Solución analítica: u(x, T) = exp(i*T) * sin(x)
        cd u_analitica = exp(cd(0, 1) * T) * sin(x_i);
        
        // Diferencia al cuadrado (norma compleja)
        suma_error_sq += pow(abs(U(i) - u_analitica), 2);
    }
    double error_global = sqrt(suma_error_sq);

    // --- EXPORTAR A .DAT ---
    ofstream outfile("resultados_p1_m1000.dat");
    if (outfile.is_open()) {
        outfile << "# x real_u imag_u abs_u" << endl;
        outfile << fixed << setprecision(8); 
        for (int i = 0; i < M + 2; ++i) {
            outfile << a + i * h << " " 
                    << U(i).real() << " " 
                    << U(i).imag() << " " 
                    << abs(U(i)) << endl;
        }
        outfile.close();
        cout << "Simulación finalizada. Datos guardados en 'resultados_p1_m1000.dat'" << endl;
    }

    cout << "------------------------------------" << endl;
    cout << "Tiempo de cálculo: " << chrono::duration<double>(end - start).count() << " s" << endl;
    cout << "Error Global (Norma L2): " << scientific << setprecision(6) << error_global << endl;
    cout << "------------------------------------" << endl;

    return 0;
}