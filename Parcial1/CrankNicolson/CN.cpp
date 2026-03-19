// ============================================================
//  Solución de la Ecuación de Difusión — Método de Crank-Nicolson
//  Traducción directa del notebook Python:
//    Crank_Nicolson_Difusion.ipynb
//
//  Problema:
//    dT/dt = alpha^2 * d²T/dx²     x in (0,l), t > 0
//
//  Condiciones de frontera:
//    T(0,t)         = 0            [Dirichlet]
//    dT/dx |_{x=l}  = 0            [Neumann — frontera aislada]
//
//  Condición inicial:
//    T(x,0) = 100 °C
//
//  Esquema CN:  A * T^{n+1} = B * T^n
//
//  Frontera Neumann: T_{N+1} = (4*T_N - T_{N-1}) / 3
//  Fila N de A: -r/3 * T_{N-1}^{n+1} + (1 + r/3) * T_N^{n+1}
//  Fila N de B:  r/3 * T_{N-1}^n    + (1 - r/3) * T_N^n
//
//  Sistema tridiagonal resuelto con Algoritmo de Thomas (O(N))
//  equivalente a la factorizacion LU del notebook (lu_factor/lu_solve).
//
//  Salida:
//    CN.dat  ->  t | T1 | T2 | ... | T10 | T_frontera
//    Consola ->  parametros, tablas comparativas, metricas
//
//  Compilar:
//    g++ -O2 -std=c++17 -o crank_nicolson crank_nicolson.cpp -lm
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
#include <algorithm>
#include <numeric>

// ── Alias de tipos ───────────────────────────────────────────
using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

// ============================================================
//  CELDA 1 — Parametros del problema
//  (equivalente al primer bloque de código del notebook)
//
//  alpha  = 0.2        Constante de difusion termica
//  L      = 1.0        Ancho de la pared [m]
//  N      = 10         Numero de divisiones espaciales
//  dx     = L / N      Paso espacial
//  T0     = 100.0      Temperatura inicial [gC]
//  dt     = 0.05       Paso de tiempo [s]  (mismo que MELIN h=0.05)
//  t_end  = 30.0       Tiempo final [s]
//  NPT    = t_end/dt   Numero de pasos de tiempo
//  r      = alpha^2*dt/dx^2   Numero de Fourier
// ============================================================
const double ALPHA = 0.2;
const double L_PAR = 1.0;
const int    N_PAR = 10;
const double DX    = L_PAR / N_PAR;
const double T0    = 100.0;
const double DT    = 0.05;
const double T_END = 30.0;
const int    NPT   = static_cast<int>(T_END / DT);   // 600 pasos

// ============================================================
//  CELDA 4 — Solucion analitica (referencia)
//
//  T(x,t) = (400/pi) * sum_{n=0}^{N_terms-1}
//              1/(2n+1) * sin( (2n+1)*pi*x / (2*L) )
//              * exp( -alpha^2 * ((2n+1)*pi/(2*L))^2 * t )
//
//  Equivalente a la funcion T_analitica(x, t, ...) del notebook.
// ============================================================
const double PI = std::acos(-1.0);

double T_analitica(double x, double t, int n_terms = 50)
{
    double suma = 0.0;
    for (int n = 0; n < n_terms; ++n) {
        double k = (2*n + 1) * PI / (2.0 * L_PAR);
        suma += (1.0 / (2*n + 1))
              * std::sin(k * x)
              * std::exp(-ALPHA*ALPHA * k*k * t);
    }
    return (400.0 / PI) * suma;
}

// ============================================================
//  CELDA 3 — Construccion de las matrices del sistema
//            tridiagonal (funcion construir_matrices)
//
//  Devuelve las diagonales de A y B:
//    aA, bA, cA  ->  sub, diagonal, super  de A (implicita)
//    aB, bB, cB  ->  sub, diagonal, super  de B (explicita)
//
//  Tamano del sistema: N+1  (indices 0 .. N)
//
//  Fila 0       Dirichlet: bA[0]=1, bB[0]=0
//  Filas 1..N-1 CN estandar
//  Fila N       Neumann:   aA[N]=-r/3, bA[N]=1+r/3
//                          aB[N]= r/3, bB[N]=1-r/3
// ============================================================
struct MatrizCN {
    Vec aA, bA, cA;   // diagonales de A (implicita)
    Vec aB, bB, cB;   // diagonales de B (explicita)
};

MatrizCN construir_matrices(int n, double r)
{
    int sz = n + 1;
    MatrizCN m;
    m.aA.assign(sz, 0.0);  m.bA.assign(sz, 0.0);  m.cA.assign(sz, 0.0);
    m.aB.assign(sz, 0.0);  m.bB.assign(sz, 0.0);  m.cB.assign(sz, 0.0);

    // ── Fila 0: Dirichlet T(0,t) = 0 ────────────────────────
    m.bA[0] = 1.0;
    m.bB[0] = 0.0;   // RHS = 0 → T_0^{n+1} = 0 siempre

    // ── Filas interiores i = 1 ... N-1 ──────────────────────
    for (int i = 1; i < n; ++i) {
        m.aA[i] = -r/2.0;   m.bA[i] = 1.0 + r;   m.cA[i] = -r/2.0;
        m.aB[i] =  r/2.0;   m.bB[i] = 1.0 - r;   m.cB[i] =  r/2.0;
    }

    // ── Fila N: Neumann (dT/dx = 0 en x = l) ────────────────
    //  Sustituyendo T_{N+1} = (4*T_N - T_{N-1})/3 en CN:
    //    Implicito:  -r/3 * T_{N-1}^{n+1} + (1 + r/3) * T_N^{n+1}
    //    Explicito:   r/3 * T_{N-1}^n     + (1 - r/3) * T_N^n
    m.aA[n] = -r/3.0;   m.bA[n] = 1.0 + r/3.0;   m.cA[n] = 0.0;
    m.aB[n] =  r/3.0;   m.bB[n] = 1.0 - r/3.0;   m.cB[n] = 0.0;

    return m;
}

// ── Imprime las primeras 5x5 de una matriz tridiagonal ───────
void imprimir_5x5(const std::string& nombre,
                  const Vec& a, const Vec& b, const Vec& c, int n)
{
    int lim = std::min(n + 1, 5);
    std::cout << nombre << " — primeras " << lim << "x" << lim << ":\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < lim; ++i) {
        std::cout << "  [";
        for (int j = 0; j < lim; ++j) {
            double val = 0.0;
            if (j == i-1) val = a[i];
            if (j == i  ) val = b[i];
            if (j == i+1) val = c[i];
            std::cout << std::setw(8) << val;
        }
        std::cout << "  ]\n";
    }
    std::cout << "\n";
}

// ============================================================
//  Algoritmo de Thomas
//  Resuelve el sistema tridiagonal A*x = d en O(N).
//  Equivalente a lu_factor / lu_solve de scipy.linalg.
//
//  a : subdiagonal   (a[0] no se usa)
//  b : diagonal principal
//  c : superdiagonal (c[sz-1] no se usa)
//  d : termino independiente (se pasa por valor, se modifica)
// ============================================================
Vec thomas(const Vec& a, const Vec& b, const Vec& c, Vec d)
{
    int sz = static_cast<int>(d.size());
    Vec c2(sz, 0.0), x(sz, 0.0);

    // Barrido hacia adelante (eliminacion)
    c2[0] = c[0] / b[0];
    d[0]  = d[0] / b[0];
    for (int i = 1; i < sz; ++i) {
        double denom = b[i] - a[i] * c2[i-1];
        c2[i] = c[i] / denom;
        d[i]  = (d[i] - a[i] * d[i-1]) / denom;
    }

    // Sustitucion hacia atras
    x[sz-1] = d[sz-1];
    for (int i = sz-2; i >= 0; --i)
        x[i] = d[i] - c2[i] * x[i+1];

    return x;
}

// ============================================================
//  Multiplicacion tridiagonal: res = M * v
//  Equivalente a:  rhs = B @ T_actual  del notebook
// ============================================================
Vec mult_tridiag(const Vec& a, const Vec& b, const Vec& c, const Vec& v)
{
    int sz = static_cast<int>(v.size());
    Vec res(sz, 0.0);
    for (int i = 0; i < sz; ++i) {
        if (i > 0)    res[i] += a[i] * v[i-1];
                      res[i] += b[i] * v[i];
        if (i < sz-1) res[i] += c[i] * v[i+1];
    }
    return res;
}

// ============================================================
//  PROGRAMA PRINCIPAL
// ============================================================
int main()
{
    double r = ALPHA*ALPHA * DT / (DX*DX);

    // ==========================================================
    //  CELDA 1 — Impresion de parametros
    // ==========================================================
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "alpha  = " << ALPHA                           << "\n";
    std::cout << "l      = " << L_PAR                          << " m\n";
    std::cout << "N      = " << N_PAR                          << "  nodos\n";
    std::cout << "dx     = " << DX                             << "\n";
    std::cout << "dt     = " << DT                             << "\n";
    std::cout << "NPT    = " << NPT                            << " pasos de tiempo\n";
    std::cout << "r      = alpha^2*dt/dx^2 = "                 << r << "\n";
    std::cout << "\n";
    std::cout << "Crank-Nicolson es incondicionalmente estable -> cualquier r es valido\n";
    std::cout << "(r = " << r << ", equivalente al MELIN: r = "
              << ALPHA*ALPHA * 0.05 / (DX*DX) << ")\n";

    // ── Grid espacial ─────────────────────────────────────────
    Vec x(N_PAR+1);
    for (int i = 0; i <= N_PAR; ++i) x[i] = i * DX;

    // ==========================================================
    //  CELDA 3 — Construccion de matrices + impresion 5x5
    // ==========================================================
    MatrizCN M = construir_matrices(N_PAR, r);

    std::cout << "\n";
    imprimir_5x5("Matriz A (implicita)", M.aA, M.bA, M.cA, N_PAR);
    imprimir_5x5("Matriz B (explicita)", M.aB, M.bB, M.cB, N_PAR);

    // ==========================================================
    //  CELDA 5 — Integracion temporal — Bucle Crank-Nicolson
    //
    //  Equivalente al notebook:
    //    lu, piv = lu_factor(A)
    //    for k in range(1, len(t_vals)):
    //        rhs        = B @ T_actual
    //        rhs[0]     = 0.0
    //        T_nuevo    = lu_solve((lu, piv), rhs)
    //        T_nuevo[0] = 0.0
    //        T_hist[:, k] = T_nuevo
    //        T_actual   = T_nuevo
    // ==========================================================

    // ── Condicion inicial ──────────────────────────────────────
    Vec T_cn(N_PAR+1, T0);
    T_cn[0] = 0.0;   // Dirichlet T(0,0) = 0

    // ── Almacenamiento completo: T_hist[nodo][paso] ───────────
    int n_pasos = NPT + 1;                          // incluye t=0
    Mat T_hist(N_PAR+1, Vec(n_pasos, 0.0));
    for (int i = 0; i <= N_PAR; ++i) T_hist[i][0] = T_cn[i];

    Vec t_vals(n_pasos);
    for (int k = 0; k < n_pasos; ++k) t_vals[k] = k * DT;

    // ── Cronometro — inicio ───────────────────────────────────
    auto t_inicio_cn = std::chrono::high_resolution_clock::now();

    // ── Bucle temporal ────────────────────────────────────────
    Vec T_actual = T_cn;
    for (int k = 1; k < n_pasos; ++k) {
        Vec rhs = mult_tridiag(M.aB, M.bB, M.cB, T_actual);
        rhs[0]  = 0.0;                               // Dirichlet
        Vec T_nuevo = thomas(M.aA, M.bA, M.cA, rhs);
        T_nuevo[0]  = 0.0;                           // reforzar Dirichlet
        for (int i = 0; i <= N_PAR; ++i) T_hist[i][k] = T_nuevo[i];
        T_actual = T_nuevo;
    }

    // ── Cronometro — fin ──────────────────────────────────────
    auto   t_fin_cn  = std::chrono::high_resolution_clock::now();
    double t_ejec_cn = std::chrono::duration<double>(t_fin_cn - t_inicio_cn).count();

    // ── Temperatura en la frontera aislada: T_{N+1} ───────────
    //  Equivalente a:  T_frontera_cn = (4*T_hist[N,:] - T_hist[N-1,:])/3
    Vec T_frontera_cn(n_pasos);
    for (int k = 0; k < n_pasos; ++k)
        T_frontera_cn[k] = (4.0*T_hist[N_PAR][k] - T_hist[N_PAR-1][k]) / 3.0;

    std::cout << "\nIntegracion completada: " << n_pasos << " pasos\n";
    std::cout << "=========================================\n";
    std::cout << std::setprecision(6);
    std::cout << "  Tiempo de ejecucion : " << t_ejec_cn << " segundos\n";
    std::cout << "  Archivo generado    : CNCPP.dat\n";
    std::cout << "=========================================\n";

    // ==========================================================
    //  CELDA 6 — Guardado del archivo CN.dat
    //
    //  Formato identico a MELINP.dat:
    //    t | T1 | T2 | ... | T10 | T_frontera
    //
    //  Equivalente a:
    //    datos_cn = np.column_stack([t_vals, T_hist[:N,:].T, T_frontera_cn])
    //    np.savetxt('CN.dat', datos_cn, fmt='%10.4f', header=header)
    // ==========================================================
    {
        std::ofstream archivo("CNCPP.dat");
        if (!archivo.is_open()) {
            std::cerr << "Error: no se pudo abrir CNCPP.dat\n";
            return 1;
        }

        // Encabezado
        archivo << std::setw(10) << "t";
        for (int i = 1; i <= N_PAR; ++i)
            archivo << std::setw(10) << ("T" + std::to_string(i));
        archivo << std::setw(12) << "T_frontera\n";

        // Datos
        archivo << std::fixed << std::setprecision(4);
        for (int k = 0; k < n_pasos; ++k) {
            archivo << std::setw(10) << t_vals[k];
            for (int i = 1; i <= N_PAR; ++i)
                archivo << std::setw(10) << T_hist[i][k];
            archivo << std::setw(12) << T_frontera_cn[k] << "\n";
        }
        archivo.close();
    }

    // ── Vista previa de las primeras 5 filas (display del notebook) ──
    std::cout << "\nPrimeras 5 filas de CNCPP.dat:\n";
    std::cout << std::setw(10) << "t";
    for (int i = 1; i <= N_PAR; ++i)
        std::cout << std::setw(10) << ("T" + std::to_string(i));
    std::cout << std::setw(12) << "T_frontera\n";
    std::cout << std::string(10 + 10*N_PAR + 12, '-') << "\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int k = 0; k < 5; ++k) {
        std::cout << std::setw(10) << t_vals[k];
        for (int i = 1; i <= N_PAR; ++i)
            std::cout << std::setw(10) << T_hist[i][k];
        std::cout << std::setw(12) << T_frontera_cn[k] << "\n";
    }

    // ==========================================================
    //  CELDA 8 — Tabla comparativa en tiempos representativos
    //  t = 0.05, 5.0, 30.0 s
    //
    //  Equivalente al notebook:
    //    for t_comp in [0.05, 5.0, 30.0]:
    //        idx   = np.argmin(np.abs(t_vals - t_comp))
    //        T_an  = T_analitica(x, t_vals[idx])
    //        T_cn_v = T_hist[:, idx]
    //        err   = np.abs(T_cn_v - T_an)
    // ==========================================================
    std::vector<double> t_comp_list = {0.05, 5.0, 30.0};

    for (double t_comp : t_comp_list) {
        // argmin( |t_vals - t_comp| )
        int idx = 0;
        double min_diff = std::abs(t_vals[0] - t_comp);
        for (int k = 1; k < n_pasos; ++k) {
            double diff = std::abs(t_vals[k] - t_comp);
            if (diff < min_diff) { min_diff = diff; idx = k; }
        }

        std::cout << "\n" << std::string(57, '=') << "\n";
        std::cout << "  t = " << std::fixed << std::setprecision(4)
                  << t_vals[idx] << " s\n";
        std::cout << std::string(57, '=') << "\n";
        std::cout << std::setw(10) << "x [m]"
                  << std::setw(14) << "CN [gC]"
                  << std::setw(16) << "Analitica [gC]"
                  << std::setw(14) << "|Error| [gC]" << "\n";
        std::cout << std::string(57, '-') << "\n";

        for (int i = 0; i <= N_PAR; ++i) {
            double xi     = x[i];
            double T_cn_v = T_hist[i][idx];
            double T_an   = T_analitica(xi, t_vals[idx]);
            double err    = std::abs(T_cn_v - T_an);
            std::cout << std::fixed << std::setprecision(4)
                      << std::setw(10) << xi
                      << std::setw(14) << T_cn_v
                      << std::setw(16) << T_an
                      << std::setw(14) << std::setprecision(6) << err << "\n";
        }
    }

    // ==========================================================
    //  CELDA 9 — Metricas globales de error
    //
    //  Calcula MAE, RMSE y error maximo sobre TODOS los nodos
    //  y TODOS los pasos de tiempo (identico al notebook).
    //
    //  Equivalente a:
    //    T_analitica_grid = zeros_like(T_hist)
    //    for k, t in enumerate(t_vals):
    //        T_analitica_grid[:, k] = T_analitica(x, t)
    //    error_abs  = abs(T_hist - T_analitica_grid)
    //    mae_global  = mean(error_abs)
    //    rmse_global = sqrt(mean((T_hist - T_analitica_grid)^2))
    //    max_error   = max(error_abs)
    //    mae_nodo    = mean(error_abs, axis=1)
    // ==========================================================
    Vec mae_nodo(N_PAR+1, 0.0);
    double mae_global  = 0.0;
    double rmse_global = 0.0;
    double max_error   = 0.0;
    long long total    = 0LL;

    for (int k = 0; k < n_pasos; ++k) {
        for (int i = 0; i <= N_PAR; ++i) {
            double T_an = T_analitica(x[i], t_vals[k]);
            double err  = std::abs(T_hist[i][k] - T_an);
            mae_nodo[i]  += err;
            mae_global   += err;
            rmse_global  += err * err;
            if (err > max_error) max_error = err;
            ++total;
        }
    }
    mae_global  /= total;
    rmse_global  = std::sqrt(rmse_global / total);
    for (int i = 0; i <= N_PAR; ++i) mae_nodo[i] /= n_pasos;

    // ── Tabla MAE por nodo ────────────────────────────────────
    std::cout << "\n\nMAE promedio por nodo (sobre todos los pasos de tiempo):\n";
    std::cout << std::string(38, '-') << "\n";
    std::cout << std::setw(10) << "x [m]"
              << std::setw(22) << "MAE promedio [gC]" << "\n";
    std::cout << std::string(38, '-') << "\n";
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i <= N_PAR; ++i)
        std::cout << std::setw(10) << x[i]
                  << std::setw(22) << mae_nodo[i] << "\n";

    std::cout << "\nError MAE global (todos los nodos y tiempos): "
              << mae_global << " gC\n";

    // ── Resumen final — cuadro de metricas ───────────────────
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════╗\n";
    std::cout << "║     RESUMEN — MÉTODO DE CRANK-NICOLSON       ║\n";
    std::cout << "╠══════════════════════════════════════════════╣\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "║  alpha (difusividad)   = " << std::setw(6) << ALPHA
              << "              ║\n";
    std::cout << "║  N (nodos espaciales)  = " << std::setw(6) << N_PAR
              << "           ║\n";
    std::cout << "║  dx                    = " << std::setw(6) << DX
              << "           ║\n";
    std::cout << "║  dt                    = " << std::setw(6) << DT
              << "          ║\n";
    std::cout << "║  r = alpha^2*dt/dx^2   = " << std::setw(6) << r
              << "          ║\n";
    std::cout << "║  Pasos de tiempo       = " << std::setw(6) << NPT
              << "           ║\n";
    std::cout << "║  Estabilidad           = incondicional      ║\n";
    std::cout << "╠══════════════════════════════════════════════╣\n";
    std::cout << std::setprecision(6);
    std::cout << "║  Tiempo de ejecucion   = " << std::setw(10) << t_ejec_cn
              << " s    ║\n";
    std::cout << "║  MAE global            = " << std::setw(10) << mae_global
              << " gC  ║\n";
    std::cout << "║  RMSE global           = " << std::setw(10) << rmse_global
              << " gC  ║\n";
    std::cout << "║  Error maximo          = " << std::setw(10) << max_error
              << " gC  ║\n";
    std::cout << "╚══════════════════════════════════════════════╝\n";

    return 0;
}