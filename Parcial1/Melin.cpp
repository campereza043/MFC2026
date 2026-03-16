// ============================================================
//  MELIN.cpp  —  Método de Líneas para la ecuación de difusión
//  Equivalente exacto del Melin.f90, con:
//    · Cronómetro de alta resolución (std::chrono)
//    · Salida MELIN_cpp.dat compatible con el notebook Python
// ============================================================
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>

// ── Parámetros del problema ──────────────────────────────────
constexpr int    N     = 10;      // divisiones espaciales
constexpr int    NPT   = 600;     // pasos de tiempo a guardar
constexpr double H     = 0.050;   // paso de integración Δt
constexpr double ALFA  = 0.2;     // difusividad térmica α
constexpr double ANCHO = 1.0;     // longitud de la pared [m]

// ── Subroutine DERIVS ────────────────────────────────────────
void derivs(const std::vector<double>& Y,
            std::vector<double>&       DYDT,
            double alfa, double dx)
{
    const double a2  = alfa * alfa;
    const double dx2 = dx * dx;

    DYDT[0] = 0.0;

    for (int i = 1; i < N - 1; ++i)
        DYDT[i] = a2 * (Y[i+1] - 2.0*Y[i] + Y[i-1]) / dx2;

    double FR_DER = (4.0*Y[N-1] - Y[N-2]) / 3.0;
    DYDT[N-1] = a2 * (FR_DER - 2.0*Y[N-1] + Y[N-2]) / dx2;
}

// ── Subroutine RK4 ───────────────────────────────────────────
void rk4(const std::vector<double>& Y,
         const std::vector<double>& DYDT,
         double t, double h,
         std::vector<double>&       YOUT,
         double alfa, double dx)
{
    const double hh = h * 0.5;
    const double h6 = h / 6.0;

    std::vector<double> YT(N), DYT(N), DYM(N);

    for (int i = 0; i < N; ++i) YT[i] = Y[i] + hh * DYDT[i];
    derivs(YT, DYT, alfa, dx);

    for (int i = 0; i < N; ++i) YT[i] = Y[i] + hh * DYT[i];
    derivs(YT, DYM, alfa, dx);

    for (int i = 0; i < N; ++i) {
        YT[i]  = Y[i] + h * DYM[i];
        DYM[i] = DYT[i] + DYM[i];
    }
    derivs(YT, DYT, alfa, dx);

    for (int i = 0; i < N; ++i)
        YOUT[i] = Y[i] + h6 * (DYDT[i] + DYT[i] + 2.0*DYM[i]);
}

// ── Escribir fila en el .dat ─────────────────────────────────
void escribir_fila(std::ofstream& archivo,
                   double t,
                   const std::vector<double>& Y)
{
    double FR_DER = (4.0*Y[N-1] - Y[N-2]) / 3.0;
    archivo << std::fixed << std::setprecision(4)
            << std::setw(10) << t;
    for (int i = 0; i < N; ++i)
        archivo << std::setw(10) << Y[i];
    archivo << std::setw(10) << FR_DER << "\n";
}

// ── Main ─────────────────────────────────────────────────────
int main()
{
    const double DX = ANCHO / N;

    std::vector<double> Y(N), DYDT(N), YOUT(N);

    for (int i = 1; i < N; ++i) Y[i] = 100.0;
    Y[0] = 0.0;

    std::ofstream dat("MELIN_cpp.dat");
    if (!dat.is_open()) {
        std::cerr << "Error: no se pudo crear MELIN_cpp.dat\n";
        return 1;
    }

    double T = 0.0;

    auto t_inicio = std::chrono::high_resolution_clock::now();

    escribir_fila(dat, T, Y);

    derivs(Y, DYDT, ALFA, DX);
    rk4(Y, DYDT, T, H, YOUT, ALFA, DX);
    T += H;
    escribir_fila(dat, T, YOUT);

    for (int j = 1; j < NPT; ++j) {
        for (int i = 0; i < N; ++i) Y[i] = YOUT[i];
        T += H;
        derivs(Y, DYDT, ALFA, DX);
        rk4(Y, DYDT, T, H, YOUT, ALFA, DX);
        escribir_fila(dat, T, YOUT);
    }

    auto t_fin = std::chrono::high_resolution_clock::now();
    double t_ejecucion =
        std::chrono::duration<double>(t_fin - t_inicio).count();

    dat.close();

    std::cout << "\n=========================================\n";
    std::cout << "  MELIN C++ — Ecuación de difusión 1D\n";
    std::cout << "-----------------------------------------\n";
    std::cout << "  α      = " << ALFA  << "\n";
    std::cout << "  l      = " << ANCHO << " m\n";
    std::cout << "  N      = " << N     << "\n";
    std::cout << "  Δx     = " << DX    << " m\n";
    std::cout << "  h      = " << H     << " s\n";
    std::cout << "  Pasos  = " << NPT   << "\n";
    std::cout << "  t_fin  = " << T     << " s\n";
    std::cout << "-----------------------------------------\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Tiempo de ejecucion: "
              << t_ejecucion << " segundos\n";
    std::cout << "  Archivo generado   : MELIN_cpp.dat\n";
    std::cout << "=========================================\n\n";

    return 0;
}