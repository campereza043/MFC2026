// ============================================================
//  Solucion de la Ecuacion de Difusion — Metodo de Crank-Nicolson
//  Traduccion literal del notebook: CN.ipynb
//
//  CELDA 1 — Parametros del problema
//  CELDA 2 — construir_matrices(N, r) : matrices densas A y B
//  CELDA 3 — Condicion inicial, lu_factor, bucle temporal,
//             T_frontera_cn, tiempo de ejecucion
//  CELDA 4 — Guardar CNPY.dat, vista previa de 5 filas
//
//  Problema:
//    dT/dt = alpha^2 * d2T/dx2     x en (0,l), t > 0
//
//  Condiciones de frontera:
//    T(0,t)        = 0             [Dirichlet]
//    dT/dx |_{x=l} = 0            [Neumann — frontera aislada]
//
//  Condicion inicial:
//    T(x,0) = 100 gC
//
//  Esquema CN — matrices DENSAS (igual que el notebook):
//    A * T^{n+1} = B * T^n
//    A factorizada UNA sola vez fuera del bucle (lu_factor),
//    y reutilizada en cada paso (lu_solve).
//    Aqui se implementa la factorizacion LU con pivoteo parcial
//    equivalente exacta a scipy.linalg.lu_factor / lu_solve.
//
//  Frontera Neumann O(Dx^2):
//    T_{N+1} = (4*T_N - T_{N-1}) / 3
//    Fila N de A: A[N][N-1]=-r/3,  A[N][N]=1+r/3
//    Fila N de B: B[N][N-1]= r/3,  B[N][N]=1-r/3
//
//  Salida:
//    CNPY.dat  ->  t | T1 | T2 | ... | T10 | T_frontera
//    Consola   ->  parametros, matrices 5x5, tiempo, 5 filas
//
//  Compilar:
//    g++ -O2 -std=c++17 -o cn_cpp CN.cpp -lm
//  Ejecutar:
//    ./cn_cpp
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
//  CELDA 2 — construir_matrices(N, r)
//
//  Construye las matrices A (implicita) y B (explicita) del
//  sistema tridiagonal de Crank-Nicolson para la ecuacion de
//  difusion 1D como MATRICES DENSAS, igual que en el notebook.
//
//  Frontera x=0  : Dirichlet T=0  -> fila 0 de identidad
//  Nodos int.    : i=1,...,N-1    -> esquema CN estandar
//  Frontera x=l  : Neumann dT/dx=0
//                  usando T_{N+1}=(4T_N-T_{N-1})/3
//
//  Equivalente exacto de:
//    A = np.zeros((size, size))
//    B = np.zeros((size, size))
//    A[0, 0] = 1.0 ; B[0, 0] = 0.0
//    for i in range(1, N):
//        A[i,i-1]=-r/2; A[i,i]=1+r; A[i,i+1]=-r/2
//        B[i,i-1]= r/2; B[i,i]=1-r; B[i,i+1]= r/2
//    A[N,N-1]=-r/3; A[N,N]=1+r/3
//    B[N,N-1]= r/3; B[N,N]=1-r/3
// ============================================================
void construir_matrices(int N, double r, Mat& A, Mat& B)
{
    int size = N + 1;   // indices 0..N

    // -- A = np.zeros((size, size))
    // -- B = np.zeros((size, size))
    A.assign(size, Vec(size, 0.0));
    B.assign(size, Vec(size, 0.0));

    // -- Fila 0: Dirichlet T(0,t) = 0
    A[0][0] = 1.0;
    B[0][0] = 0.0;   // RHS siempre 0 -> T_0^{n+1} = 0

    // -- Filas interiores i = 1 ... N-1
    for (int i = 1; i < N; ++i) {
        A[i][i-1] = -r / 2.0;
        A[i][i  ] =  1.0 + r;
        A[i][i+1] = -r / 2.0;

        B[i][i-1] =  r / 2.0;
        B[i][i  ] =  1.0 - r;
        B[i][i+1] =  r / 2.0;
    }

    // -- Fila N: Neumann (dT/dx = 0 en x = l)
    // Sustituyendo T_{N+1} = (4*T_N - T_{N-1})/3 en el esquema CN:
    // Lado implicito:
    //   -r/2*T_{N-1}^{n+1} + (1+r)*T_N^{n+1} - r/2*(4T_N-T_{N-1})/3
    //   = -r/3*T_{N-1}^{n+1} + (1+r/3)*T_N^{n+1}
    A[N][N-1] = -r / 3.0;
    A[N][N  ] =  1.0 + r / 3.0;

    // Lado explicito (analogo):
    B[N][N-1] =  r / 3.0;
    B[N][N  ] =  1.0 - r / 3.0;
}

// ============================================================
//  Imprime las primeras 5x5 de una matriz densa.
//  Equivalente a:  print(np.round(A[:5, :5], 4))
// ============================================================
void imprimir_5x5(const std::string& nombre, const Mat& M)
{
    int lim = std::min((int)M.size(), 5);
    std::cout << nombre << "\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int i = 0; i < lim; ++i) {
        std::cout << "  [";
        for (int j = 0; j < lim; ++j) {
            // imprimir -0.0 como 0.0 (estetica igual a numpy)
            double v = M[i][j];
            if (v == 0.0) v = 0.0;
            std::cout << std::setw(7) << v;
            if (j < lim - 1) std::cout << " ";
        }
        std::cout << " ]\n";
    }
    std::cout << "\n";
}

// ============================================================
//  lu_factor equivalente a scipy.linalg.lu_factor
//  Factorizacion LU con pivoteo parcial en la propia matriz.
//  Devuelve la matriz LU (in-place en A_lu) y el vector piv.
//
//  Convencion identica a scipy:
//    - La triangular inferior L (sin diagonal) y la superior U
//      quedan solapadas en A_lu.
//    - piv[i] = indice de la fila con la que se intercambio la
//      fila i durante la eliminacion.
// ============================================================
void lu_factor(Mat A_lu, Mat& LU, Vec& piv_d, std::vector<int>& piv)
{
    int n = (int)A_lu.size();
    LU  = A_lu;
    piv.resize(n);

    for (int k = 0; k < n; ++k) {
        // -- Buscar pivote maximo en la columna k (pivoteo parcial)
        int max_idx = k;
        double max_val = std::abs(LU[k][k]);
        for (int i = k+1; i < n; ++i) {
            if (std::abs(LU[i][k]) > max_val) {
                max_val = std::abs(LU[i][k]);
                max_idx = i;
            }
        }
        piv[k] = max_idx;

        // -- Intercambiar filas k y max_idx
        if (max_idx != k)
            std::swap(LU[k], LU[max_idx]);

        // -- Eliminacion gaussiana
        if (LU[k][k] != 0.0) {
            for (int i = k+1; i < n; ++i) {
                LU[i][k] /= LU[k][k];
                for (int j = k+1; j < n; ++j)
                    LU[i][j] -= LU[i][k] * LU[k][j];
            }
        }
    }
    (void)piv_d;   // no se usa, solo para mantener firma similar
}

// ============================================================
//  lu_solve equivalente a scipy.linalg.lu_solve((lu, piv), rhs)
//  Resuelve A*x = b dado la factorizacion LU con pivoteo.
// ============================================================
Vec lu_solve(const Mat& LU, const std::vector<int>& piv, Vec b)
{
    int n = (int)b.size();

    // -- Aplicar permutaciones de fila al RHS
    for (int k = 0; k < n; ++k)
        if (piv[k] != k)
            std::swap(b[k], b[piv[k]]);

    // -- Sustitucion hacia adelante: L * y = b
    for (int i = 1; i < n; ++i)
        for (int j = 0; j < i; ++j)
            b[i] -= LU[i][j] * b[j];

    // -- Sustitucion hacia atras: U * x = y
    for (int i = n-1; i >= 0; --i) {
        for (int j = i+1; j < n; ++j)
            b[i] -= LU[i][j] * b[j];
        b[i] /= LU[i][i];
    }

    return b;
}

// ============================================================
//  Multiplicacion matriz densa por vector: res = M * v
//  Equivalente a:  rhs = B @ T_actual
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
    std::cout << std::setprecision(1);
    std::cout << "dx     = " << dx << "\n";
    std::cout << std::setprecision(2);
    std::cout << "dt     = " << dt << "\n";
    std::cout << "NPT    = " << NPT << " pasos de tiempo\n";
    std::cout << std::setprecision(4);
    std::cout << "r      = alpha^2*dt/dx^2 = " << r << "\n";
    std::cout << "\n";
    std::cout << "Crank-Nicolson es incondicionalmente estable -> cualquier r es valido\n";
    std::cout << "(r = " << r << ", equivalente al MELIN: r = "
              << alpha*alpha * 0.05 / (dx*dx) << ")\n";

    // ==========================================================
    //  CELDA 2 — construir_matrices(N, r)
    //
    //  A, B = construir_matrices(N, r)
    //  print("Matriz A (implicita) — primeras 5x5:")
    //  print(np.round(A[:5, :5], 4))
    //  print()
    //  print("Matriz B (explicita) — primeras 5x5:")
    //  print(np.round(B[:5, :5], 4))
    // ==========================================================
    Mat A, B;
    construir_matrices(N, r, A, B);

    std::cout << "\nMatriz A (implicita) -- primeras 5x5:\n";
    imprimir_5x5("", A);

    std::cout << "Matriz B (explicita) -- primeras 5x5:\n";
    imprimir_5x5("", B);

    // ==========================================================
    //  CELDA 3 — Condicion inicial, lu_factor, bucle temporal
    //
    //  T_cn        = np.full(N+1, T0)
    //  T_cn[0]     = 0.0
    //  t_vals      = np.arange(0, t_end + dt, dt)
    //  T_hist      = np.zeros((N+1, len(t_vals)))
    //  T_hist[:,0] = T_cn.copy()
    //  lu, piv     = lu_factor(A)          <- UNA sola vez
    //  t_inicio_cn = time.perf_counter()
    //  for k in range(1, len(t_vals)):
    //      rhs        = B @ T_actual
    //      rhs[0]     = 0.0
    //      T_nuevo    = lu_solve((lu, piv), rhs)
    //      T_nuevo[0] = 0.0
    //      T_hist[:,k]= T_nuevo
    //      T_actual   = T_nuevo
    //  t_ejec_cn = time.perf_counter() - t_inicio_cn
    //  T_frontera_cn = (4*T_hist[N,:] - T_hist[N-1,:])/3.0
    // ==========================================================

    // -- Condicion inicial
    // T_cn = np.full(N+1, T0)
    Vec T_cn(N+1, T0);
    T_cn[0] = 0.0;   // Dirichlet T(0,0) = 0

    // -- Almacenamiento de la solucion
    // t_vals = np.arange(0, t_end + dt, dt)   -> NPT+1 valores
    int n_tvals = NPT + 1;
    Vec t_vals(n_tvals);
    for (int k = 0; k < n_tvals; ++k) t_vals[k] = k * dt;

    // T_hist = np.zeros((N+1, len(t_vals)))
    // Almacenamos como T_hist[nodo][paso] para acceso por filas
    Mat T_hist(N+1, Vec(n_tvals, 0.0));

    // T_hist[:, 0] = T_cn.copy()
    for (int i = 0; i <= N; ++i) T_hist[i][0] = T_cn[i];

    // -- Factorizacion LU de A (solo una vez, fuera del bucle)
    // lu, piv = lu_factor(A)
    Mat LU;
    Vec piv_dummy;
    std::vector<int> piv;
    lu_factor(A, LU, piv_dummy, piv);

    // -- Cronometro: t_inicio_cn = time.perf_counter()
    auto t_inicio_cn = std::chrono::high_resolution_clock::now();

    // -- Bucle temporal
    // T_actual = T_cn.copy()
    Vec T_actual = T_cn;

    for (int k = 1; k < n_tvals; ++k) {
        // rhs = B @ T_actual
        Vec rhs = mat_vec(B, T_actual);

        // rhs[0] = 0.0   (Dirichlet: T_0 siempre 0)
        rhs[0] = 0.0;

        // T_nuevo = lu_solve((lu, piv), rhs)
        Vec T_nuevo = lu_solve(LU, piv, rhs);

        // T_nuevo[0] = 0.0   (reforzar Dirichlet)
        T_nuevo[0] = 0.0;

        // T_hist[:, k] = T_nuevo
        for (int i = 0; i <= N; ++i) T_hist[i][k] = T_nuevo[i];

        // T_actual = T_nuevo
        T_actual = T_nuevo;
    }

    // -- t_ejec_cn = time.perf_counter() - t_inicio_cn
    auto t_fin_cn   = std::chrono::high_resolution_clock::now();
    double t_ejec_cn = std::chrono::duration<double>(t_fin_cn - t_inicio_cn).count();

    // -- Temperatura en la frontera aislada (x=l)
    // T_frontera_cn = (4*T_hist[N, :] - T_hist[N-1, :]) / 3.0
    Vec T_frontera_cn(n_tvals);
    for (int k = 0; k < n_tvals; ++k)
        T_frontera_cn[k] = (4.0*T_hist[N][k] - T_hist[N-1][k]) / 3.0;

    // -- Impresion del resultado de la celda 3
    std::cout << "Integracion completada: " << n_tvals << " pasos\n";
    std::cout << "=========================================\n";
    std::cout << std::setprecision(6);
    std::cout << "  Tiempo de ejecucion : " << t_ejec_cn << " segundos\n";
    std::cout << "  Archivo generado    : CNCPP.dat\n";
    std::cout << "=========================================\n";

    // ==========================================================
    //  CELDA 4 — Guardar CNPY.dat + vista previa de 5 filas
    //
    //  datos_cn = np.column_stack([
    //      t_vals,
    //      T_hist[:N, :].T,    # nodos 1..N (sin nodo 0, siempre 0)
    //      T_frontera_cn
    //  ])
    //  header = "t" + "".join([f"      T{i+1}" for i in range(N)])
    //           + "  T_frontera"
    //  np.savetxt('CNPY.dat', datos_cn, fmt='%10.4f',
    //             header=header, comments='')
    //
    //  print("Primeras 5 filas de CNPY.dat:")
    //  display(df_preview.round(4))
    // ==========================================================

    // -- Abrir CNPY.dat y escribir encabezado
    std::ofstream archivo("CNCPP.dat");
    if (!archivo.is_open()) {
        std::cerr << "Error: no se pudo crear CNCPP.dat\n";
        return 1;
    }

    // header = "t" + "".join([f"      T{i+1}" for i in range(N)]) + "  T_frontera"
    // El Python guarda T_hist[:N, :].T que son los indices Python 0..N-1,
    // es decir los nodos T0 (Dirichlet=0) hasta T_{N-1}, mas T_frontera.
    // Por eso la primera columna de datos es T0 = 0.0 siempre.
    archivo << std::setw(10) << "t";
    for (int i = 0; i < N; ++i)
        archivo << std::setw(10) << ("T" + std::to_string(i+1));
    archivo << std::setw(12) << "T_frontera" << "\n";

    // datos: t | T0..T_{N-1} | T_frontera   (fmt='%10.4f')
    // T_hist[:N, :].T  ->  nodos indice Python 0..N-1
    //                   = nodos C++ T_hist[0]..T_hist[N-1]
    archivo << std::fixed << std::setprecision(4);
    for (int k = 0; k < n_tvals; ++k) {
        archivo << std::setw(10) << t_vals[k];
        for (int i = 0; i < N; ++i)          // nodos Python 0..N-1
            archivo << std::setw(10) << T_hist[i][k];
        archivo << std::setw(12) << T_frontera_cn[k] << "\n";
    }
    archivo.close();

    // -- Vista previa de las primeras 5 filas
    // print("Primeras 5 filas de CNPY.dat:")
    // display(df_preview.round(4))
    std::cout << "\nPrimeras 5 filas de CNCPP.dat:\n";
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
        std::cout << std::setw(12) << T_frontera_cn[k] << "\n";
    }

    return 0;
}