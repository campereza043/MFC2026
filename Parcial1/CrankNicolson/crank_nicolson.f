C ============================================================
C  Solucion de la Ecuacion de Difusion - Metodo de Crank-Nicolson
C  Traduccion directa del notebook Python:
C    Crank_Nicolson_Difusion.ipynb
C
C  Problema:
C    dT/dt = alpha^2 * d2T/dx2     x en (0,l), t > 0
C
C  Condiciones de frontera:
C    T(0,t)        = 0             [Dirichlet]
C    dT/dx|_{x=l} = 0             [Neumann - frontera aislada]
C
C  Condicion inicial:
C    T(x,0) = 100 gC
C
C  Esquema CN:  A * T^{n+1} = B * T^n
C
C  Frontera Neumann O(Dx^2):
C    T_{N+1} = (4*T_N - T_{N-1}) / 3
C  Fila N de A: -r/3*T_{N-1}^{n+1} + (1+r/3)*T_N^{n+1}
C  Fila N de B:  r/3*T_{N-1}^n     + (1-r/3)*T_N^n
C
C  Sistema tridiagonal resuelto con Algoritmo de Thomas (O(N))
C  equivalente a lu_factor/lu_solve de scipy.linalg.
C
C  Salida:
C    CN.dat  ->  t | T1 | T2 | ... | T10 | T_frontera
C    Consola ->  parametros, matrices 5x5, tablas, metricas
C
C  Compilar:
C    gfortran -O2 -o crank_nicolson_f crank_nicolson.f
C  Ejecutar:
C    ./crank_nicolson_f
C ============================================================

      PROGRAM CRANK_NICOLSON
      IMPLICIT REAL*8 (A-H, O-Z)

C ============================================================
C  CELDA 1 - Parametros del problema
C  (equivalente al primer bloque de codigo del notebook)
C
C  ALPHA  = 0.2        Constante de difusion termica
C  ANCHO  = 1.0        Ancho de la pared [m]
C  N      = 10         Numero de divisiones espaciales
C  DX     = ANCHO/N    Paso espacial
C  T0     = 100.0      Temperatura inicial [gC]
C  DT     = 0.05       Paso de tiempo [s] (mismo que MELIN h=0.05)
C  T_END  = 30.0       Tiempo final [s]
C  NPT    = T_END/DT   Numero de pasos de tiempo
C  R      = ALPHA^2*DT/DX^2   Numero de Fourier
C ============================================================

      PARAMETER (N=10, NPT=600, NMAX=N+1, NPMAX=NPT+1)

C     -- Variables del problema
      REAL*8 ALPHA, ANCHO, DX, T0, DT, T_END, R
      REAL*8 PI

C     -- Diagonales de A (implicita) y B (explicita)
      REAL*8 AA(NMAX), BA(NMAX), CA(NMAX)
      REAL*8 AB(NMAX), BB(NMAX), CB(NMAX)

C     -- Vectores de trabajo
      REAL*8 X(NMAX)
      REAL*8 T_ACTUAL(NMAX), T_NUEVO(NMAX)
      REAL*8 RHS(NMAX)

C     -- Almacenamiento completo T_HIST(nodo, paso)
      REAL*8 T_HIST(NMAX, NPMAX)
      REAL*8 T_VALS(NPMAX)
      REAL*8 T_FRONT(NPMAX)

C     -- Metricas de error
      REAL*8 MAE_NODO(NMAX)
      REAL*8 MAE_GLOBAL, RMSE_GLOBAL, MAX_ERROR

C     -- Temporales de calculo
      REAL*8 T_AN, ERR, SUMA_MAE, SUMA_RMSE
      REAL*8 T_EJEC_CN, T_COMP

C     -- Variables para el cronometro (CPU time)
      REAL*8 CPU_INI, CPU_FIN

C     -- Indices
      INTEGER I, K, PASO, IDX, TOTAL, NPASOS

C     -- Parametros fisicos
      ALPHA = 0.2D0
      ANCHO = 1.0D0
      DX    = ANCHO / DBLE(N)
      T0    = 100.0D0
      DT    = 0.05D0
      T_END = 30.0D0
      PI    = DACOS(-1.0D0)
      R     = ALPHA*ALPHA * DT / (DX*DX)
      NPASOS = NPT + 1

C     -- Impresion de parametros (equivalente al print del notebook)
      WRITE(*,'(A,F6.1)')   'alpha  = ', ALPHA
      WRITE(*,'(A,F6.1,A)') 'l      = ', ANCHO, ' m'
      WRITE(*,'(A,I3,A)')   'N      = ', N,     '  nodos'
      WRITE(*,'(A,F6.4)')   'dx     = ', DX
      WRITE(*,'(A,F6.4)')   'dt     = ', DT
      WRITE(*,'(A,I4,A)')   'NPT    = ', NPT,   ' pasos de tiempo'
      WRITE(*,'(A,F8.4)')   'r      = alpha^2*dt/dx^2 = ', R
      WRITE(*,*)
      WRITE(*,'(A)') 'Crank-Nicolson es incondicionalmente estable'
      WRITE(*,'(A,F6.4,A,F6.4,A)')
     &  '(r = ', R, ', equivalente al MELIN: r = ',
     &  ALPHA*ALPHA*0.05D0/(DX*DX), ')'

C     -- Grid espacial: x_i = i * DX   i = 0,...,N
      DO I = 1, NMAX
         X(I) = DBLE(I-1) * DX
      END DO

C ============================================================
C  CELDA 3 - Construccion de las matrices del sistema
C            tridiagonal (funcion construir_matrices)
C
C  Indices Fortran: 1..NMAX  corresponden a  0..N del Python
C
C  Fila 1   (i=0 Python): Dirichlet BA(1)=1, BB(1)=0
C  Filas 2..N (i=1..N-1): CN estandar
C  Fila NMAX (i=N Python): Neumann
C    AA(NMAX)=-r/3, BA(NMAX)=1+r/3
C    AB(NMAX)= r/3, BB(NMAX)=1-r/3
C ============================================================

C     -- Inicializar a cero
      DO I = 1, NMAX
         AA(I) = 0.0D0
         BA(I) = 0.0D0
         CA(I) = 0.0D0
         AB(I) = 0.0D0
         BB(I) = 0.0D0
         CB(I) = 0.0D0
      END DO

C     -- Fila 1: Dirichlet T(0,t) = 0
      BA(1) = 1.0D0
      BB(1) = 0.0D0

C     -- Filas interiores i = 2 ... N  (Python i=1..N-1)
      DO I = 2, N
         AA(I) = -R/2.0D0
         BA(I) =  1.0D0 + R
         CA(I) = -R/2.0D0
         AB(I) =  R/2.0D0
         BB(I) =  1.0D0 - R
         CB(I) =  R/2.0D0
      END DO

C     -- Fila NMAX: Neumann (dT/dx = 0 en x = l)
C       Sustituyendo T_{N+1} = (4*T_N - T_{N-1})/3 en CN:
C         Implicito: -r/3 * T_{N-1}^{n+1} + (1+r/3) * T_N^{n+1}
C         Explicito:  r/3 * T_{N-1}^n     + (1-r/3) * T_N^n
      AA(NMAX) = -R/3.0D0
      BA(NMAX) =  1.0D0 + R/3.0D0
      CA(NMAX) =  0.0D0
      AB(NMAX) =  R/3.0D0
      BB(NMAX) =  1.0D0 - R/3.0D0
      CB(NMAX) =  0.0D0

C     -- Imprimir primeras 5x5 de A y B
C       (equivalente al print de numpy del notebook)
      WRITE(*,*)
      WRITE(*,'(A)') 'Matriz A (implicita) - primeras 5x5:'
      CALL IMPR5X5(AA, BA, CA, NMAX)

      WRITE(*,'(A)') 'Matriz B (explicita) - primeras 5x5:'
      CALL IMPR5X5(AB, BB, CB, NMAX)

C ============================================================
C  CELDA 5 - Integracion temporal - Bucle Crank-Nicolson
C
C  Equivalente al notebook:
C    lu, piv = lu_factor(A)
C    for k in range(1, len(t_vals)):
C        rhs        = B @ T_actual
C        rhs[0]     = 0.0
C        T_nuevo    = lu_solve((lu, piv), rhs)
C        T_nuevo[0] = 0.0
C        T_hist[:, k] = T_nuevo
C        T_actual   = T_nuevo
C ============================================================

C     -- Condicion inicial
C       T_cn = np.full(N+1, T0);  T_cn[0] = 0.0
      DO I = 1, NMAX
         T_ACTUAL(I) = T0
      END DO
      T_ACTUAL(1) = 0.0D0

C     -- Almacenar t_vals y T_hist en paso 0
      DO K = 1, NPASOS
         T_VALS(K) = DBLE(K-1) * DT
      END DO
      DO I = 1, NMAX
         T_HIST(I, 1) = T_ACTUAL(I)
      END DO

C     -- Cronometro: inicio (CPU_TIME equivale a time.perf_counter)
      CALL CPU_TIME(CPU_INI)

C     -- Bucle temporal: k = 2 .. NPASOS  (Python k = 1 .. len-1)
      DO PASO = 2, NPASOS

C        -- RHS = B * T_ACTUAL  (mult_tridiag del C++ / B @ T_actual Python)
         CALL MULT_TRIDIAG(AB, BB, CB, T_ACTUAL, RHS, NMAX)

C        -- Dirichlet en RHS: rhs[0] = 0
         RHS(1) = 0.0D0

C        -- Resolver A * T_NUEVO = RHS con Algoritmo de Thomas
C           (equivalente a lu_solve del notebook)
         CALL THOMAS(AA, BA, CA, RHS, T_NUEVO, NMAX)

C        -- Reforzar Dirichlet en solucion
         T_NUEVO(1) = 0.0D0

C        -- Guardar en T_HIST
         DO I = 1, NMAX
            T_HIST(I, PASO) = T_NUEVO(I)
         END DO

C        -- Avanzar en el tiempo
         DO I = 1, NMAX
            T_ACTUAL(I) = T_NUEVO(I)
         END DO

      END DO

C     -- Cronometro: fin
      CALL CPU_TIME(CPU_FIN)
      T_EJEC_CN = CPU_FIN - CPU_INI

C     -- Temperatura frontera aislada T_{N+1}
C       T_frontera_cn = (4*T_hist[N,:] - T_hist[N-1,:])/3
      DO K = 1, NPASOS
         T_FRONT(K) = (4.0D0*T_HIST(NMAX,K) - T_HIST(NMAX-1,K))/3.0D0
      END DO

      WRITE(*,*)
      WRITE(*,'(A,I4,A)') 'Integracion completada: ', NPASOS, ' pasos'
      WRITE(*,'(A)') '========================================='
      WRITE(*,'(A,F12.6,A)') '  Tiempo de ejecucion : ',
     &                        T_EJEC_CN, ' segundos'
      WRITE(*,'(A)') '  Archivo generado    : CN.dat'
      WRITE(*,'(A)') '========================================='

C ============================================================
C  CELDA 6 - Guardado del archivo CN.dat
C
C  Formato identico a MELINP.dat:
C    t | T1 | T2 | ... | T10 | T_frontera
C
C  Equivalente a:
C    datos_cn = np.column_stack([t_vals, T_hist[:N,:].T, T_frontera_cn])
C    np.savetxt('CN.dat', datos_cn, fmt='%10.4f', header=header)
C ============================================================

      OPEN(UNIT=2, FILE='CN.dat', STATUS='UNKNOWN')

C     -- Encabezado
      WRITE(2,'(A10,10A10,A12)')
     &  't','T1','T2','T3','T4','T5',
     &  'T6','T7','T8','T9','T10','T_frontera'

C     -- Datos: t | T1..T10 | T_frontera
      DO K = 1, NPASOS
         WRITE(2,'(11F10.4, F12.4)')
     &     T_VALS(K),
     &     (T_HIST(I,K), I=2,NMAX),
     &     T_FRONT(K)
      END DO

      CLOSE(2)

C     -- Vista previa de las primeras 5 filas (display del notebook)
      WRITE(*,*)
      WRITE(*,'(A)') 'Primeras 5 filas de CN.dat:'
      WRITE(*,'(A10,10A10,A12)')
     &  't','T1','T2','T3','T4','T5',
     &  'T6','T7','T8','T9','T10','T_frontera'
      WRITE(*,'(A)') REPEAT('-', 122)
      DO K = 1, 5
         WRITE(*,'(11F10.4, F12.4)')
     &     T_VALS(K),
     &     (T_HIST(I,K), I=2,NMAX),
     &     T_FRONT(K)
      END DO

C ============================================================
C  CELDA 8 - Tabla comparativa en tiempos representativos
C  t = 0.05, 5.0, 30.0 s
C
C  Equivalente al notebook:
C    for t_comp in [0.05, 5.0, 30.0]:
C        idx   = np.argmin(np.abs(t_vals - t_comp))
C        T_an  = T_analitica(x, t_vals[idx])
C        T_cn_v = T_hist[:, idx]
C        err   = np.abs(T_cn_v - T_an)
C ============================================================

      REAL*8 T_COMP_LIST(3)
      T_COMP_LIST(1) = 0.05D0
      T_COMP_LIST(2) = 5.0D0
      T_COMP_LIST(3) = 30.0D0

      DO ITCOMP = 1, 3
         T_COMP = T_COMP_LIST(ITCOMP)

C        -- argmin( |t_vals - t_comp| )
         IDX = 1
         DIFFMIN = DABS(T_VALS(1) - T_COMP)
         DO K = 2, NPASOS
            DIFF = DABS(T_VALS(K) - T_COMP)
            IF (DIFF .LT. DIFFMIN) THEN
               DIFFMIN = DIFF
               IDX = K
            END IF
         END DO

         WRITE(*,*)
         WRITE(*,'(A)') REPEAT('=', 57)
         WRITE(*,'(A,F8.4,A)') '  t = ', T_VALS(IDX), ' s'
         WRITE(*,'(A)') REPEAT('=', 57)
         WRITE(*,'(A10,A14,A16,A14)')
     &     'x [m]', 'CN [gC]', 'Analitica [gC]', '|Error| [gC]'
         WRITE(*,'(A)') REPEAT('-', 57)

         DO I = 1, NMAX
            T_AN = T_ANALITICA(X(I), T_VALS(IDX), ALPHA, ANCHO, PI)
            ERR  = DABS(T_HIST(I,IDX) - T_AN)
            WRITE(*,'(F10.4, F14.4, F16.4, F14.6)')
     &        X(I), T_HIST(I,IDX), T_AN, ERR
         END DO

      END DO

C ============================================================
C  CELDA 9 - Metricas globales de error
C
C  Calcula MAE, RMSE y error maximo sobre TODOS los nodos
C  y TODOS los pasos de tiempo (identico al notebook).
C
C  Equivalente a:
C    T_analitica_grid[:,k] = T_analitica(x, t_vals[k])
C    error_abs  = abs(T_hist - T_analitica_grid)
C    mae_global  = mean(error_abs)
C    rmse_global = sqrt(mean((T_hist - T_analitica_grid)^2))
C    max_error   = max(error_abs)
C    mae_nodo    = mean(error_abs, axis=1)
C ============================================================

      DO I = 1, NMAX
         MAE_NODO(I) = 0.0D0
      END DO
      MAE_GLOBAL  = 0.0D0
      RMSE_GLOBAL = 0.0D0
      MAX_ERROR   = 0.0D0
      TOTAL       = 0

      DO K = 1, NPASOS
         DO I = 1, NMAX
            T_AN = T_ANALITICA(X(I), T_VALS(K), ALPHA, ANCHO, PI)
            ERR  = DABS(T_HIST(I,K) - T_AN)
            MAE_NODO(I)  = MAE_NODO(I)  + ERR
            MAE_GLOBAL   = MAE_GLOBAL   + ERR
            RMSE_GLOBAL  = RMSE_GLOBAL  + ERR*ERR
            IF (ERR .GT. MAX_ERROR) MAX_ERROR = ERR
            TOTAL = TOTAL + 1
         END DO
      END DO

      MAE_GLOBAL  = MAE_GLOBAL  / DBLE(TOTAL)
      RMSE_GLOBAL = DSQRT(RMSE_GLOBAL / DBLE(TOTAL))
      DO I = 1, NMAX
         MAE_NODO(I) = MAE_NODO(I) / DBLE(NPASOS)
      END DO

C     -- Tabla MAE por nodo
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,'(A)') 'MAE promedio por nodo (sobre todos los pasos):'
      WRITE(*,'(A)') REPEAT('-', 38)
      WRITE(*,'(A10,A22)') 'x [m]', 'MAE promedio [gC]'
      WRITE(*,'(A)') REPEAT('-', 38)
      DO I = 1, NMAX
         WRITE(*,'(F10.6, F22.6)') X(I), MAE_NODO(I)
      END DO

      WRITE(*,*)
      WRITE(*,'(A,F12.6,A)') 'Error MAE global (todos nodos y tiempos):',
     &  MAE_GLOBAL, ' gC'

C     -- Cuadro de metricas final (identico al notebook)
      WRITE(*,*)
      WRITE(*,'(A)') '╔══════════════════════════════════════════════╗'
      WRITE(*,'(A)') '║     RESUMEN - METODO DE CRANK-NICOLSON       ║'
      WRITE(*,'(A)') '╠══════════════════════════════════════════════╣'
      WRITE(*,'(A,F6.1,A)')
     &  '║  alpha (difusividad)   = ', ALPHA, '              ║'
      WRITE(*,'(A,I6,A)')
     &  '║  N (nodos espaciales)  = ', N,     '           ║'
      WRITE(*,'(A,F6.4,A)')
     &  '║  dx                    = ', DX,    '           ║'
      WRITE(*,'(A,F6.4,A)')
     &  '║  dt                    = ', DT,    '          ║'
      WRITE(*,'(A,F6.4,A)')
     &  '║  r = alpha^2*dt/dx^2   = ', R,     '          ║'
      WRITE(*,'(A,I6,A)')
     &  '║  Pasos de tiempo       = ', NPT,   '           ║'
      WRITE(*,'(A)') '║  Estabilidad           = incondicional      ║'
      WRITE(*,'(A)') '╠══════════════════════════════════════════════╣'
      WRITE(*,'(A,F12.6,A)')
     &  '║  Tiempo de ejecucion   = ', T_EJEC_CN,  ' s    ║'
      WRITE(*,'(A,F12.6,A)')
     &  '║  MAE global            = ', MAE_GLOBAL,  ' gC  ║'
      WRITE(*,'(A,F12.6,A)')
     &  '║  RMSE global           = ', RMSE_GLOBAL, ' gC  ║'
      WRITE(*,'(A,F12.6,A)')
     &  '║  Error maximo          = ', MAX_ERROR,   ' gC  ║'
      WRITE(*,'(A)') '╚══════════════════════════════════════════════╝'

      STOP
      END

C ============================================================
C  SUBRUTINA IMPR5X5
C  Imprime las primeras 5x5 de una matriz tridiagonal.
C  Equivalente al print(np.round(A[:5,:5], 4)) del notebook.
C ============================================================
      SUBROUTINE IMPR5X5(A, B, C, NMAX)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER NMAX, LIM, I, J
      REAL*8  A(NMAX), B(NMAX), C(NMAX), VAL

      LIM = MIN(NMAX, 5)
      DO I = 1, LIM
         WRITE(*,'(2X,A)', ADVANCE='NO') '['
         DO J = 1, LIM
            VAL = 0.0D0
            IF (J .EQ. I-1) VAL = A(I)
            IF (J .EQ. I  ) VAL = B(I)
            IF (J .EQ. I+1) VAL = C(I)
            WRITE(*,'(F8.4)', ADVANCE='NO') VAL
         END DO
         WRITE(*,'(A)') '  ]'
      END DO
      WRITE(*,*)
      RETURN
      END

C ============================================================
C  SUBRUTINA MULT_TRIDIAG
C  Calcula RES = M * V  para matriz tridiagonal M = (A,B,C).
C  Equivalente a:  rhs = B @ T_actual   del notebook Python.
C ============================================================
      SUBROUTINE MULT_TRIDIAG(A, B, C, V, RES, NMAX)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER NMAX, I
      REAL*8  A(NMAX), B(NMAX), C(NMAX), V(NMAX), RES(NMAX)

      DO I = 1, NMAX
         RES(I) = B(I) * V(I)
         IF (I .GT. 1)    RES(I) = RES(I) + A(I) * V(I-1)
         IF (I .LT. NMAX) RES(I) = RES(I) + C(I) * V(I+1)
      END DO
      RETURN
      END

C ============================================================
C  SUBRUTINA THOMAS
C  Resuelve el sistema tridiagonal A*X = D en O(N).
C  Equivalente a lu_factor / lu_solve de scipy.linalg.
C
C  A   : subdiagonal   (A(1) no se usa)
C  B   : diagonal principal
C  C   : superdiagonal (C(NMAX) no se usa)
C  D   : termino independiente (se pasa por valor via DCOPY)
C  X   : solucion
C ============================================================
      SUBROUTINE THOMAS(A, B, C, D, X, NMAX)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER NMAX, I
      REAL*8  A(NMAX), B(NMAX), C(NMAX), D(NMAX), X(NMAX)
      REAL*8  C2(NMAX), D2(NMAX), DENOM

C     -- Copiar D en D2 para no modificar el original
      DO I = 1, NMAX
         C2(I) = 0.0D0
         D2(I) = D(I)
      END DO

C     -- Barrido hacia adelante (eliminacion)
      C2(1) = C(1) / B(1)
      D2(1) = D2(1) / B(1)
      DO I = 2, NMAX
         DENOM = B(I) - A(I) * C2(I-1)
         C2(I) = C(I) / DENOM
         D2(I) = (D2(I) - A(I) * D2(I-1)) / DENOM
      END DO

C     -- Sustitucion hacia atras
      X(NMAX) = D2(NMAX)
      DO I = NMAX-1, 1, -1
         X(I) = D2(I) - C2(I) * X(I+1)
      END DO

      RETURN
      END

C ============================================================
C  FUNCION T_ANALITICA
C  Solucion analitica por separacion de variables (50 terminos).
C  Equivalente a la funcion T_analitica(x, t, ...) del notebook:
C
C  T(x,t) = (400/pi) * sum_{n=0}^{49}
C              1/(2n+1) * sin((2n+1)*pi*x/(2*L))
C              * exp(-alpha^2 * ((2n+1)*pi/(2*L))^2 * t)
C ============================================================
      REAL*8 FUNCTION T_ANALITICA(X, T, ALPHA, L, PI)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER N_TERMS, N
      PARAMETER (N_TERMS = 50)
      REAL*8  X, T, ALPHA, L, PI, K_N, SUMA

      SUMA = 0.0D0
      DO N = 0, N_TERMS-1
         K_N  = DBLE(2*N + 1) * PI / (2.0D0 * L)
         SUMA = SUMA + (1.0D0 / DBLE(2*N + 1))
     &              * DSIN(K_N * X)
     &              * DEXP(-ALPHA*ALPHA * K_N*K_N * T)
      END DO
      T_ANALITICA = (400.0D0 / PI) * SUMA

      RETURN
      END
