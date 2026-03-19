// ============================================================
//  Solucion de la Ecuacion de Difusion
//  Metodo de Diferencias Finitas Explicito
//  Traduccion literal del notebook: EXP.ipynb
//
//  CELDA 1 - Parametros del problema
//  CELDA 2 - construir_matriz_explicita(N, r) : matriz densa B
//  CELDA 3 - Condicion inicial, bucle temporal T^{n+1}=B*T^n,
//             T_frontera, tiempo de ejecucion
//  CELDA 4 - Guardar EXP.dat, vista previa de 5 filas
//
//  Problema:
//    dT/dt = alpha^2 * d2T/dx2     x en (0,l), t > 0
//
//  Condiciones de frontera:
//    T(0,t)        = 0             [Dirichlet]
//    dT/dx |_{x=l} = 0            [Neumann - frontera aislada]
//
//  Condicion inicial:
//    T(x,0) = 100 gC
//
//  Esquema EXPLICITO:
//    T^{n+1} = B * T^n      (solo un producto matriz-vector)
//    No hay sistema lineal que resolver.
//
//    Nodos interiores:  T_i^{n+1} = r*T_{i-1} + (1-2r)*T_i + r*T_{i+1}
//    Frontera Neumann:  T_N^{n+1} = (2r/3)*T_{N-1} + (1-2r/3)*T_N
//                       (usando T_{N+1} = (4T_N - T_{N-1})/3)
//
//  Condicion de estabilidad: r = alpha^2*dt/dx^2 <= 0.5
//
//  Salida:
//    EXP.dat  ->  t | T1 | T2 | ... | T10 | T_frontera
//    Consola  ->  parametros, matriz 5x5, tiempo, 5 filas
//
//  Compilar:
//    g++ -O2 -std=c++17 -o exp_cpp EXP.cpp -lm
//  Ejecutar:
//    ./exp_cpp
// ============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <string>
#include <algorithm>

// ── Alias de tipos (matrices densas, igual que numpy) ────────
using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

// ============================================================
//  CELDA 2 — construir_matriz_explicita(N, r)
//
//  Construye la matriz explícita B como MATRIZ DENSA.
//  El paso temporal es simplemente:  T^{n+1} = B * T^n
//
//  Frontera x=0  : Dirichlet T=0  -> fila 0 de ceros
//  Nodos int.    : i=1,...,N-1    -> esquema explicito estandar
//                  B[i,i-1]=r, B[i,i]=1-2r, B[i,i+1]=r
//  Frontera x=l  : Neumann dT/dx=0
//                  usando T_{N+1}=(4T_N-T_{N-1})/3:
//                  B[N,N-1]=2r/3, B[N,N]=1-2r/3
//
//  Equivalente exacto de:
//    B = np.zeros((size, size))
//    B[0, 0] = 0.0
//    for i in range(1, N):
//        B[i,i-1]=r; B[i,i]=1-2r; B[i,i+1]=r
//    B[N,N-1]=2r/3; B[N,N]=1-2r/3
// ============================================================
void construir_matriz_explicita(int N, double r, Mat& B)
{
    int size = N + 1;   // indices 0..N

    // -- B = np.zeros((size, size))
    B.assign(size, Vec(size, 0.0));

    // -- Fila 0: Dirichlet T(0,t) = 0
    B[0][0] = 0.0;   // T_0^{n+1} = 0 siempre

    // -- Filas interiores i = 1 ... N-1
    for (int i = 1; i < N; ++i) {
        B[i][i-1] =  r;
        B[i][i  ] =  1.0 - 2.0*r;
        B[i][i+1] =  r;
    }

    // -- Fila N: Neumann (dT/dx = 0 en x = l)
    // Sustituyendo T_{N+1} = (4*T_N - T_{N-1})/3 en el esquema:
    //   T_N^{n+1} = r*T_{N-1} + (1-2r)*T_N + r*(4T_N-T_{N-1})/3
    //             = (r - r/3)*T_{N-1} + (1-2r + 4r/3)*T_N
    //             = (2r/3)*T_{N-1}   + (1 - 2r/3)*T_N
    B[N][N-1] =  2.0*r / 3.0;
    B[N][N  ] =  1.0 - 2.0*r / 3.0;
}

// ============================================================
//  Imprime las primeras 5x5 de una matriz densa.
//  Equivalente a:  print(np.round(B[:5, :5], 4))
// ============================================================
void imprimir_5x5(const Mat& M)
{
    int lim = std::min((int)M.size(), 5);
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < lim; ++i) {
        std::cout << "  [";
        for (int j = 0; j < lim; ++j) {
            std::cout << std::setw(7) << M[i][j];
            if (j < lim - 1) std::cout << " ";
        }
        std::cout << " ]\n";
    }
    std::cout << "\n";
}

// ============================================================
//  Multiplicacion matriz densa por vector: res = M * v
//  Equivalente a:  T_nuevo = B @ T_actual
// ============================================================
Vec mat_vec(const Mat& M, const Vec& v)
{
    int n = (int)v.size();
    Vec res(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            res[i] += M[i][j] * v[j];
    return res;
}

// ============================================================
//  PROGRAMA PRINCIPAL
// ============================================================
int main()
{
    // ==========================================================
    //  CELDA 1 — Parametros del problema
    //  (equivalente al primer bloque de codigo del notebook)
    // ==========================================================

    // -- Parametros del problema (identicos a MELIN)
    const double alpha = 0.2;       // Constante de difusion termica
    const double L     = 1.0;       // Ancho de la pared [m]
    const int    N     = 10;        // Numero de divisiones espaciales
    const double dx    = L / N;     // Paso espacial
    const double T0    = 100.0;     // Temperatura inicial [gC]

    // -- Parametros temporales
    const double dt    = 0.05;              // Paso de tiempo [s] (mismo que MELIN h=0.05)
    const double t_end = 30.0;             // Tiempo final [s]
    const int    NPT   = (int)(t_end/dt);  // Numero de pasos de tiempo

    // -- Numero de Fourier (parametro de estabilidad)
    const double r = alpha*alpha * dt / (dx*dx);

    // -- Grid espacial: x = np.linspace(0, L, N+1)
    Vec x(N+1);
    for (int i = 0; i <= N; ++i) x[i] = i * dx;

    // -- Impresion de parametros
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "alpha  = " << alpha << "\n";
    std::cout << "l      = " << L << " m\n";
    std::cout << "N      = " << N << "  nodos\n";
    std::cout << "dx     = " << dx << "\n";
    std::cout << std::setprecision(2);
    std::cout << "dt     = " << dt << "\n";
    std::cout << "NPT    = " << NPT << " pasos de tiempo\n";
    std::cout << std::setprecision(4);
    std::cout << "r      = alpha^2*dt/dx^2 = " << r << "\n";
    std::cout << "\n";
    if (r <= 0.5)
        std::cout << "Criterio de estabilidad: r = " << r << " <= 0.5 -> ESTABLE\n";
    else
        std::cout << "Criterio de estabilidad: r = " << r << " >  0.5 -> INESTABLE\n";

    // ==========================================================
    //  CELDA 2 — construir_matriz_explicita(N, r)
    //
    //  B = construir_matriz_explicita(N, r)
    //  print("Matriz B (explicita) — primeras 5x5:")
    //  print(np.round(B[:5, :5], 4))
    // ==========================================================
    Mat B;
    construir_matriz_explicita(N, r, B);

    std::cout << "\nMatriz B (explicita) -- primeras 5x5:\n";
    imprimir_5x5(B);

    // ==========================================================
    //  CELDA 3 — Condicion inicial, bucle temporal
    //
    //  T_exp       = np.full(N+1, T0)
    //  T_exp[0]    = 0.0
    //  t_vals      = np.arange(0, t_end + dt, dt)
    //  T_hist      = np.zeros((N+1, len(t_vals)))
    //  T_hist[:,0] = T_exp.copy()
    //  t_inicio    = time.perf_counter()
    //  for k in range(1, len(t_vals)):
    //      T_nuevo    = B @ T_actual    <- avance explicito
    //      T_nuevo[0] = 0.0
    //      T_hist[:,k]= T_nuevo
    //      T_actual   = T_nuevo
    //  t_ejec = time.perf_counter() - t_inicio
    //  T_frontera = (4*T_hist[N,:] - T_hist[N-1,:])/3.0
    // ==========================================================

    // -- Condicion inicial
    // T_exp = np.full(N+1, T0)
    Vec T_exp(N+1, T0);
    T_exp[0] = 0.0;   // Dirichlet T(0,0) = 0

    // -- Almacenamiento de la solucion
    // t_vals = np.arange(0, t_end + dt, dt)   -> NPT+1 valores
    int n_tvals = NPT + 1;
    Vec t_vals(n_tvals);
    for (int k = 0; k < n_tvals; ++k) t_vals[k] = k * dt;

    // T_hist = np.zeros((N+1, len(t_vals)))
    Mat T_hist(N+1, Vec(n_tvals, 0.0));

    // T_hist[:, 0] = T_exp.copy()
    for (int i = 0; i <= N; ++i) T_hist[i][0] = T_exp[i];

    // -- Cronometro: t_inicio = time.perf_counter()
    auto t_inicio = std::chrono::high_resolution_clock::now();

    // -- Bucle temporal: T^{n+1} = B * T^n
    // T_actual = T_exp.copy()
    Vec T_actual = T_exp;

    for (int k = 1; k < n_tvals; ++k) {
        // T_nuevo = B @ T_actual   (avance explicito)
        Vec T_nuevo = mat_vec(B, T_actual);

        // T_nuevo[0] = 0.0   (reforzar Dirichlet)
        T_nuevo[0] = 0.0;

        // T_hist[:, k] = T_nuevo
        for (int i = 0; i <= N; ++i) T_hist[i][k] = T_nuevo[i];

        // T_actual = T_nuevo
        T_actual = T_nuevo;
    }

    // -- t_ejec = time.perf_counter() - t_inicio
    auto t_fin  = std::chrono::high_resolution_clock::now();
    double t_ejec = std::chrono::duration<double>(t_fin - t_inicio).count();

    // -- Temperatura en la frontera aislada (x=l)
    // T_frontera = (4*T_hist[N, :] - T_hist[N-1, :]) / 3.0
    Vec T_frontera(n_tvals);
    for (int k = 0; k < n_tvals; ++k)
        T_frontera[k] = (4.0*T_hist[N][k] - T_hist[N-1][k]) / 3.0;

    // -- Impresion del resultado de la celda 3
    std::cout << "Integracion completada: " << n_tvals << " pasos\n";
    std::cout << "=========================================\n";
    std::cout << std::setprecision(6);
    std::cout << "  Tiempo de ejecucion : " << t_ejec << " segundos\n";
    std::cout << "  Archivo generado    : FTCSCPP.dat\n";
    std::cout << "=========================================\n";

    // ==========================================================
    //  CELDA 4 — Guardar EXP.dat + vista previa de 5 filas
    //
    //  datos_exp = np.column_stack([
    //      t_vals,
    //      T_hist[:N, :].T,   # nodos Python 0..N-1 (T0=0 incluido)
    //      T_frontera
    //  ])
    //  header = "t      T1 ... T10  T_frontera"
    //  np.savetxt('EXP.dat', datos_exp, fmt='%10.4f', ...)
    //
    //  print("Primeras 5 filas de EXP.dat:")
    //  display(df_preview.round(4))
    //
    //  Nota: T_hist[:N, :].T = nodos Python 0..N-1
    //                        = C++ T_hist[0]..T_hist[N-1]
    //    Primera columna de datos = T_hist[0][k] = 0 (Dirichlet)
    // ==========================================================

    // -- Abrir EXP.dat y escribir encabezado
    std::ofstream archivo("FTCSCPP.dat");
    if (!archivo.is_open()) {
        std::cerr << "Error: no se pudo crear FTCSCPP.dat\n";
        return 1;
    }

    // header = "t" + "".join([f"      T{i+1}" for i in range(N)]) + "  T_frontera"
    archivo << std::setw(10) << "t";
    for (int i = 0; i < N; ++i)
        archivo << std::setw(10) << ("T" + std::to_string(i+1));
    archivo << std::setw(12) << "T_frontera" << "\n";

    // datos: t | T_hist[0..N-1][k] | T_frontera[k]   (fmt='%10.4f')
    archivo << std::fixed << std::setprecision(4);
    for (int k = 0; k < n_tvals; ++k) {
        archivo << std::setw(10) << t_vals[k];
        for (int i = 0; i < N; ++i)          // nodos Python 0..N-1
            archivo << std::setw(10) << T_hist[i][k];
        archivo << std::setw(12) << T_frontera[k] << "\n";
    }
    archivo.close();

    // -- Vista previa de las primeras 5 filas
    std::cout << "\nPrimeras 5 filas de FTCSCPP.dat:\n";
    std::cout << std::setw(10) << "t";
    for (int i = 0; i < N; ++i)
        std::cout << std::setw(10) << ("T" + std::to_string(i+1));
    std::cout << std::setw(12) << "T_frontera" << "\n";
    std::cout << std::string(10 + 10*N + 12, '-') << "\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int k = 0; k < 5; ++k) {
        std::cout << std::setw(10) << t_vals[k];
        for (int i = 0; i < N; ++i)
            std::cout << std::setw(10) << T_hist[i][k];
        std::cout << std::setw(12) << T_frontera[k] << "\n";
    }

    return 0;
}