// ============================================================
//  SOLUCIÓN DE LA ECUACIÓN DE DIFUSIÓN
//  MÉTODO DE CRANK-NICOLSON — C++
// ============================================================
//
//  Problema:
//    ∂T/∂t = α² ∂²T/∂x²    x ∈ (0,l), t > 0
//
//  Condiciones de frontera:
//    T(0,t)            = 0       (Dirichlet)
//    ∂T/∂x |_{x=l}    = 0       (Neumann — frontera aislada)
//
//  Condición inicial:
//    T(x,0) = 100 °C
//
//  Frontera Neumann tratada con diferencia backward O(Δx²):
//    T_{N+1} = (4*T_N - T_{N-1}) / 3
//
//  Sistema tridiagonal en cada paso de tiempo:
//    A * T^{n+1} = B * T^n
//
//  Compilar:
//    g++ -O2 -o crank_nicolson crank_nicolson.cpp
//  Ejecutar:
//    ./crank_nicolson
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <string>

// ── Constantes del problema ──────────────────────────────────
const double ALPHA  = 0.2;     // Difusividad térmica
const double L      = 1.0;     // Ancho de la pared [m]
const int    N      = 10;      // Número de divisiones espaciales
const double DX     = L / N;   // Paso espacial
const double T0     = 100.0;   // Temperatura inicial [°C]
const double DT     = 0.05;    // Paso de tiempo [s]
const double T_END  = 30.0;    // Tiempo final [s]
const int    NPT    = static_cast<int>(T_END / DT);  // Pasos de tiempo

// ── Tipo alias para matrices y vectores ──────────────────────
using Vec  = std::vector<double>;
using Mat  = std::vector<Vec>;

// ============================================================
//  thomas_algorithm
//  Resuelve un sistema tridiagonal  A*x = d
//  mediante el algoritmo de Thomas (eliminación gaussiana
//  optimizada para matrices tridiagonales — O(N)).
//
//  Entradas:
//    a : subdiagonal  (a[0] no se usa)
//    b : diagonal principal
//    c : superdiagonal (c[N-1] no se usa)
//    d : término independiente
//  Salida:
//    x : solución del sistema
// ============================================================
Vec thomas_algorithm(const Vec& a, const Vec& b, const Vec& c, Vec d)
{
    int n = static_cast<int>(d.size());
    Vec x(n, 0.0);
    Vec c_prima(n, 0.0);
    Vec d_prima(n, 0.0);

    // ── Barrido hacia adelante ────────────────────────────
    c_prima[0] = c[0] / b[0];
    d_prima[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m   = a[i] / (b[i] - a[i] * c_prima[i-1]);
        c_prima[i] = c[i] / (b[i] - a[i] * c_prima[i-1]);
        d_prima[i] = (d[i] - a[i] * d_prima[i-1])
                     / (b[i] - a[i] * c_prima[i-1]);
    }

    // ── Sustitución hacia atrás ───────────────────────────
    x[n-1] = d_prima[n-1];
    for (int i = n-2; i >= 0; --i)
        x[i] = d_prima[i] - c_prima[i] * x[i+1];

    return x;
}

// ============================================================
//  T_analitica
//  Solución analítica por separación de variables:
//
//  T(x,t) = (400/π) Σ_{n=0}^{N_terms-1}
//              1/(2n+1) * sin((2n+1)π x/(2l))
//              * exp(-α² ((2n+1)π/(2l))² t)
// ============================================================
double T_analitica(double x, double t, int n_terms = 50)
{
    double suma = 0.0;
    for (int n = 0; n < n_terms; ++n) {
        double k = (2*n + 1) * M_PI / (2.0 * L);
        suma += (1.0 / (2*n + 1))
              * std::sin(k * x)
              * std::exp(-ALPHA*ALPHA * k*k * t);
    }
    return (400.0 / M_PI) * suma;
}

// ============================================================
//  construir_diagonales
//  Arma las diagonales de las matrices A (implícita) y B
//  (explícita) del esquema de Crank-Nicolson.
//
//  Tamaño del sistema: N+1 nodos (índices 0 .. N)
//
//  Fila 0   → Dirichlet: T_0 = 0
//  Filas 1..N-1 → CN estándar
//  Fila N   → Neumann:  T_{N+1} = (4T_N - T_{N-1})/3
//               ⟹ -r/3 * T_{N-1}^{n+1} + (1+r/3)*T_N^{n+1} = …
// ============================================================
struct Tridiagonal {
    Vec a;   // subdiagonal
    Vec b;   // diagonal
    Vec c;   // superdiagonal
};

struct SistemasCN {
    Tridiagonal A;   // lado implícito
    Tridiagonal B;   // lado explícito
};

SistemasCN construir_diagonales(int n, double r)
{
    int sz = n + 1;
    SistemasCN sys;

    sys.A.a.assign(sz, 0.0);
    sys.A.b.assign(sz, 0.0);
    sys.A.c.assign(sz, 0.0);

    sys.B.a.assign(sz, 0.0);
    sys.B.b.assign(sz, 0.0);
    sys.B.c.assign(sz, 0.0);

    // ── Fila 0: Dirichlet ─────────────────────────────────
    sys.A.b[0] = 1.0;
    sys.B.b[0] = 0.0;   // RHS = 0 → T_0^{n+1} = 0

    // ── Filas interiores ──────────────────────────────────
    for (int i = 1; i < n; ++i) {
        sys.A.a[i] = -r/2.0;
        sys.A.b[i] =  1.0 + r;
        sys.A.c[i] = -r/2.0;

        sys.B.a[i] =  r/2.0;
        sys.B.b[i] =  1.0 - r;
        sys.B.c[i] =  r/2.0;
    }

    // ── Fila N: Neumann ───────────────────────────────────
    sys.A.a[n] = -r/3.0;
    sys.A.b[n] =  1.0 + r/3.0;
    sys.A.c[n] =  0.0;

    sys.B.a[n] =  r/3.0;
    sys.B.b[n] =  1.0 - r/3.0;
    sys.B.c[n] =  0.0;

    return sys;
}

// ============================================================
//  multiplicar_tridiagonal
//  Multiplica una matriz tridiagonal (a,b,c) por un vector v:
//  resultado[i] = a[i]*v[i-1] + b[i]*v[i] + c[i]*v[i+1]
// ============================================================
Vec multiplicar_tridiagonal(const Tridiagonal& M, const Vec& v)
{
    int n = static_cast<int>(v.size());
    Vec res(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (i > 0)   res[i] += M.a[i] * v[i-1];
                     res[i] += M.b[i] * v[i];
        if (i < n-1) res[i] += M.c[i] * v[i+1];
    }
    return res;
}

// ============================================================
//  PROGRAMA PRINCIPAL
// ============================================================
int main()
{
    // ── Número de Fourier ────────────────────────────────────
    double r = ALPHA*ALPHA * DT / (DX*DX);

    // ── Encabezado en consola ────────────────────────────────
    std::cout << "=============================================" << std::endl;
    std::cout << "  ECUACIÓN DE DIFUSIÓN — CRANK-NICOLSON C++  " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  α     = " << ALPHA  << std::endl;
    std::cout << "  l     = " << L      << " m"  << std::endl;
    std::cout << "  N     = " << N      << " nodos" << std::endl;
    std::cout << "  Δx    = " << DX     << std::endl;
    std::cout << "  Δt    = " << DT     << std::endl;
    std::cout << "  NPT   = " << NPT    << " pasos de tiempo" << std::endl;
    std::cout << "  r     = α²Δt/Δx² = " << r << std::endl;
    std::cout << "  Estabilidad: incondicional (CN) ✔" << std::endl;
    std::cout << "=============================================" << std::endl;

    // ── Grid espacial ────────────────────────────────────────
    Vec x(N+1);
    for (int i = 0; i <= N; ++i)
        x[i] = i * DX;

    // ── Condición inicial ────────────────────────────────────
    Vec T_actual(N+1, T0);
    T_actual[0] = 0.0;   // Dirichlet T(0,0)=0

    // ── Construir matrices ───────────────────────────────────
    SistemasCN sys = construir_diagonales(N, r);

    // ── Abrir archivo de salida CN.dat ───────────────────────
    std::ofstream archivo("CNCPP.dat");
    if (!archivo.is_open()) {
        std::cerr << "Error: no se pudo abrir CNCPP.dat" << std::endl;
        return 1;
    }

    // Encabezado del archivo
    archivo << std::setw(10) << "t";
    for (int i = 1; i <= N; ++i)
        archivo << std::setw(10) << ("T" + std::to_string(i));
    archivo << std::setw(12) << "T_frontera" << "\n";

    // ── Escribir condición inicial (t=0) ─────────────────────
    double T_frontera = (4.0*T_actual[N] - T_actual[N-1]) / 3.0;
    archivo << std::setw(10) << std::fixed << std::setprecision(4) << 0.0;
    for (int i = 1; i <= N; ++i)
        archivo << std::setw(10) << T_actual[i];
    archivo << std::setw(12) << T_frontera << "\n";

    // ── Inicio del cronómetro ────────────────────────────────
    auto t_inicio = std::chrono::high_resolution_clock::now();

    // ── Bucle temporal ───────────────────────────────────────
    double t = 0.0;
    for (int paso = 1; paso <= NPT; ++paso) {
        t += DT;

        // RHS = B * T_actual
        Vec rhs = multiplicar_tridiagonal(sys.B, T_actual);

        // Forzar Dirichlet en RHS
        rhs[0] = 0.0;

        // Resolver A * T_nuevo = rhs con algoritmo de Thomas
        Vec T_nuevo = thomas_algorithm(sys.A.a, sys.A.b, sys.A.c, rhs);

        // Forzar Dirichlet en la solución
        T_nuevo[0] = 0.0;

        // Temperatura virtual en frontera aislada
        T_frontera = (4.0*T_nuevo[N] - T_nuevo[N-1]) / 3.0;

        // Escribir en archivo
        archivo << std::setw(10) << t;
        for (int i = 1; i <= N; ++i)
            archivo << std::setw(10) << T_nuevo[i];
        archivo << std::setw(12) << T_frontera << "\n";

        // Avanzar en el tiempo
        T_actual = T_nuevo;
    }

    // ── Fin del cronómetro ───────────────────────────────────
    auto t_fin   = std::chrono::high_resolution_clock::now();
    double t_ejec = std::chrono::duration<double>(t_fin - t_inicio).count();

    archivo.close();

    // ── Resultado final en consola ───────────────────────────
    std::cout << std::endl;
    std::cout << "  Integración completada: " << NPT << " pasos" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Tiempo de ejecucion : " << t_ejec << " segundos" << std::endl;
    std::cout << "  Archivo generado    : CNCPP.dat"  << std::endl;
    std::cout << "=============================================" << std::endl;

    // ── Comparación con solución analítica en t=5 y t=30 ────
    std::cout << std::endl;
    std::cout << "  Comparación CN vs Analítica — t = 5.0 s" << std::endl;
    std::cout << "  " << std::string(52, '-') << std::endl;
    std::cout << std::setw(8)  << "x [m]"
              << std::setw(12) << "CN [°C]"
              << std::setw(14) << "Analít. [°C]"
              << std::setw(12) << "|Error|" << std::endl;
    std::cout << "  " << std::string(52, '-') << std::endl;

    // Para t=5, recalcular desde cero (usamos la solución guardada)
    // Leemos los datos desde el archivo CN.dat para t=5.0
    // En su lugar mostramos la analítica evaluada en x[i]
    // y el valor final de T_actual (que es t=30)
    // → Recalculamos T(x,5) corriendo solo los pasos necesarios
    Vec T_t5(N+1, T0);
    T_t5[0] = 0.0;
    int pasos_5 = static_cast<int>(5.0 / DT);
    for (int paso = 0; paso < pasos_5; ++paso) {
        Vec rhs = multiplicar_tridiagonal(sys.B, T_t5);
        rhs[0]  = 0.0;
        T_t5    = thomas_algorithm(sys.A.a, sys.A.b, sys.A.c, rhs);
        T_t5[0] = 0.0;
    }

    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i <= N; ++i) {
        double xi   = x[i];
        double T_an = T_analitica(xi, 5.0);
        double err  = std::abs(T_t5[i] - T_an);
        std::cout << std::setw(8)  << xi
                  << std::setw(12) << T_t5[i]
                  << std::setw(14) << T_an
                  << std::setw(12) << err << std::endl;
    }

    std::cout << std::endl;
    std::cout << "  Comparación CN vs Analítica — t = 30.0 s" << std::endl;
    std::cout << "  " << std::string(52, '-') << std::endl;
    std::cout << std::setw(8)  << "x [m]"
              << std::setw(12) << "CN [°C]"
              << std::setw(14) << "Analít. [°C]"
              << std::setw(12) << "|Error|" << std::endl;
    std::cout << "  " << std::string(52, '-') << std::endl;

    for (int i = 0; i <= N; ++i) {
        double xi   = x[i];
        double T_an = T_analitica(xi, 30.0);
        double err  = std::abs(T_actual[i] - T_an);
        std::cout << std::setw(8)  << xi
                  << std::setw(12) << T_actual[i]
                  << std::setw(14) << T_an
                  << std::setw(12) << err << std::endl;
    }

    // ── Métricas globales de error (sobre todo el dominio en t=30) ──
    double mae = 0.0, rmse = 0.0, max_err = 0.0;
    for (int i = 0; i <= N; ++i) {
        double T_an  = T_analitica(x[i], 30.0);
        double err   = std::abs(T_actual[i] - T_an);
        mae    += err;
        rmse   += err * err;
        if (err > max_err) max_err = err;
    }
    mae  /= (N + 1);
    rmse  = std::sqrt(rmse / (N + 1));

    std::cout << std::endl;
    std::cout << "╔══════════════════════════════════════════════╗" << std::endl;
    std::cout << "║     RESUMEN — MÉTODO DE CRANK-NICOLSON       ║" << std::endl;
    std::cout << "╠══════════════════════════════════════════════╣" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "║  Tiempo de ejecución  = " << t_ejec  << " s    ║" << std::endl;
    std::cout << "║  MAE   (en t=30)      = " << mae     << " °C  ║" << std::endl;
    std::cout << "║  RMSE  (en t=30)      = " << rmse    << " °C  ║" << std::endl;
    std::cout << "║  Error máximo (t=30)  = " << max_err << " °C  ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════════╝" << std::endl;

    return 0;
}