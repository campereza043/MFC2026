#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>

int main() {
    // ── Parámetros ────────────────────────────────────────────
    const int N = 10;
    const int NPT = 600;
    const double H = 0.050;
    const double ALFA = 0.2;
    const double ANCHO = 1.0;

    // ── Variables ─────────────────────────────────────────────
    // Usamos tamaño N+1 para mantener los índices de 1 a N como en Fortran
    std::vector<double> Y(N + 1), Y_NEW(N + 1);
    std::vector<double> A(N + 1), B(N + 1), C(N + 1), D(N + 1); // sistema tridiagonal

    double T, DX, R;
    double FR_DER;

    DX = ANCHO / N;
    T = 0.0;
    R = (ALFA * ALFA * H) / (2.0 * DX * DX); // número de Fourier/2

    // Abrir archivo de salida con el nuevo nombre
    std::ofstream outFile("CNCPP.dat");
    if (!outFile.is_open()) {
        std::cerr << "Error al abrir el archivo CNCPP.dat" << std::endl;
        return 1;
    }

    // ── Condición inicial ──────────────────────────────────────
    for (int i = 2; i <= N; ++i) {
        Y[i] = 100.0;
    }
    Y[1] = 0.0; // frontera fría fija

    // ── Cronómetro ────────────────────────────────────────────
    auto t_inicio = std::chrono::high_resolution_clock::now();

    // ── Escribir t=0 ──────────────────────────────────────────
    FR_DER = (4.0 * Y[N] - Y[N - 1]) / 3.0;
    
    outFile << std::fixed << std::setprecision(4);
    outFile << std::setw(10) << T;
    for (int i = 1; i <= N; ++i) {
        outFile << std::setw(10) << Y[i];
    }
    outFile << std::setw(10) << FR_DER << "\n";

    // ── Ciclo temporal ────────────────────────────────────────
    for (int j = 1; j <= NPT; ++j) {
        
        // Construir sistema tridiagonal A·Y_NEW = D
        // Frontera x=0: Dirichlet T1=0 siempre
        A[1] = 0.0;
        B[1] = 1.0;
        C[1] = 0.0;
        D[1] = 0.0;

        // Puntos interiores i=2..N-1
        for (int i = 2; i <= N - 1; ++i) {
            A[i] = -R;
            B[i] = 1.0 + 2.0 * R;
            C[i] = -R;
            D[i] = R * Y[i - 1] + (1.0 - 2.0 * R) * Y[i] + R * Y[i + 1];
        }

        // Frontera x=l: Neumann dT/dx=0
        A[N] = -R / 3.0;
        B[N] = 1.0 + (2.0 * R) / 3.0;
        C[N] = 0.0;
        D[N] = (R / 3.0) * Y[N - 1] + (1.0 - (2.0 * R) / 3.0) * Y[N];

        // ── Eliminación de Thomas (TDMA) ──────────────────────
        // Se calcula un multiplicador temporal para optimizar la división
        for (int i = 2; i <= N; ++i) {
            double m = A[i] / B[i - 1];
            B[i] = B[i] - m * C[i - 1];
            D[i] = D[i] - m * D[i - 1];
        }

        // Sustitución regresiva
        Y_NEW[N] = D[N] / B[N];
        for (int i = N - 1; i >= 1; --i) {
            Y_NEW[i] = (D[i] - C[i] * Y_NEW[i + 1]) / B[i];
        }

        // Actualización temporal
        T = T + H;
        for (int i = 1; i <= N; ++i) {
            Y[i] = Y_NEW[i];
        }

        FR_DER = (4.0 * Y[N] - Y[N - 1]) / 3.0;
        
        // Escribir paso temporal
        outFile << std::setw(10) << T;
        for (int i = 1; i <= N; ++i) {
            outFile << std::setw(10) << Y[i];
        }
        outFile << std::setw(10) << FR_DER << "\n";
    }

    // ── Fin del Cronómetro ────────────────────────────────────
    auto t_fin = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_ejecucion = t_fin - t_inicio;

    outFile.close();

    // ── Prints de salida editados ─────────────────────────────
    std::cout << "=================================\n";
    std::cout << "  MELIN — Crank-Nicolson C++\n";
    std::cout << "---------------------------------\n";
    std::cout << "  R (Fourier/2) = " << std::fixed << std::setprecision(4) << R << "\n";
    std::cout << "  Estabilidad   : INCONDICIONAL\n";
    std::cout << "  Tiempo        : " << std::fixed << std::setprecision(6) << t_ejecucion.count() << " s\n";
    std::cout << "  Archivo       : CNCPP.dat\n";
    std::cout << "=================================\n";

    return 0;
}