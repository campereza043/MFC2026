!============================================================
!  Solucion de la Ecuacion de Difusion - Metodo de Crank-Nicolson
!  Traduccion literal del notebook: CN.ipynb
!  (via CN.cpp — mismo orden, misma logica, mismas funciones)
!
!  CELDA 1 - Parametros del problema
!  CELDA 2 - construir_matrices(N, r) : matrices DENSAS A y B
!  CELDA 3 - Condicion inicial, lu_factor, bucle temporal,
!             T_frontera_cn, tiempo de ejecucion
!  CELDA 4 - Guardar CNF.dat, vista previa de 5 filas
!
!  Problema:
!    dT/dt = alpha^2 * d2T/dx2     x en (0,l), t > 0
!
!  Condiciones de frontera:
!    T(0,t)        = 0             [Dirichlet]
!    dT/dx |_{x=l} = 0            [Neumann - frontera aislada]
!
!  Condicion inicial:
!    T(x,0) = 100 gC
!
!  Esquema CN - matrices DENSAS (igual que el notebook):
!    A * T^{n+1} = B * T^n
!    A factorizada UNA sola vez fuera del bucle (lu_factor),
!    y reutilizada en cada paso (lu_solve).
!    Factorizacion LU con pivoteo parcial equivalente exacta
!    a scipy.linalg.lu_factor / lu_solve.
!
!  Frontera Neumann O(Dx^2):
!    T_{N+1} = (4*T_N - T_{N-1}) / 3
!    Fila N+1 de A: A(N+1,N)=-r/3,  A(N+1,N+1)=1+r/3
!    Fila N+1 de B: B(N+1,N)= r/3,  B(N+1,N+1)=1-r/3
!
!  Salida:
!    CNF.dat  ->  t | T1 | T2 | ... | T10 | T_frontera
!    Consola   ->  parametros, matrices 5x5, tiempo, 5 filas
!
!  Compilar:
!    gfortran -O2 -o cn_f90 CN.f90
!  Ejecutar:
!    ./cn_f90
! ============================================================


! ============================================================
!  MODULO: params
!  Centraliza precision y dimensiones (equivalente a los
!  #define / const del C++).
! ============================================================
MODULE params
  IMPLICIT NONE
  INTEGER,  PARAMETER :: dp   = KIND(1.0D0)   ! precision doble
  INTEGER,  PARAMETER :: N    = 10             ! divisiones espaciales
  INTEGER,  PARAMETER :: NPT  = 600            ! pasos de tiempo
  INTEGER,  PARAMETER :: SZ   = N + 1          ! tamano del sistema (0..N)
  INTEGER,  PARAMETER :: NTOT = NPT + 1        ! tiempos totales  (0..NPT)
END MODULE params


! ============================================================
!  MODULO: rutinas
!  Contiene todas las subrutinas en el mismo orden que el C++:
!    - construir_matrices   (CELDA 2)
!    - imprimir_5x5         (equivalente a print numpy)
!    - lu_factor_sub        (equivalente a scipy lu_factor)
!    - lu_solve_sub         (equivalente a scipy lu_solve)
!    - mat_vec              (equivalente a B @ T_actual)
! ============================================================
MODULE rutinas
  USE params
  IMPLICIT NONE

CONTAINS

  ! ----------------------------------------------------------
  !  SUBRUTINA: construir_matrices
  !  Construye A (implicita) y B (explicita) como matrices
  !  DENSAS de tamano SZ x SZ. Equivalente exacto de:
  !
  !    A = np.zeros((size, size))
  !    B = np.zeros((size, size))
  !    A[0,0]=1.0 ; B[0,0]=0.0
  !    for i in range(1, N):
  !        A[i,i-1]=-r/2; A[i,i]=1+r; A[i,i+1]=-r/2
  !        B[i,i-1]= r/2; B[i,i]=1-r; B[i,i+1]= r/2
  !    A[N,N-1]=-r/3; A[N,N]=1+r/3
  !    B[N,N-1]= r/3; B[N,N]=1-r/3
  !
  !  Nota de indices: Fortran 1..SZ <-> Python/C++ 0..N
  !    Fila 1   (Python 0)  : Dirichlet
  !    Filas 2..N (Python 1..N-1) : CN estandar
  !    Fila SZ  (Python N)  : Neumann
  ! ----------------------------------------------------------
  SUBROUTINE construir_matrices(r, A, B)
    REAL(dp), INTENT(IN)  :: r
    REAL(dp), INTENT(OUT) :: A(SZ, SZ), B(SZ, SZ)
    INTEGER :: i

    ! -- A = np.zeros((size, size)) ; B = np.zeros((size, size))
    A = 0.0_dp
    B = 0.0_dp

    ! -- Fila 1 (Python fila 0): Dirichlet T(0,t) = 0
    A(1, 1) = 1.0_dp
    B(1, 1) = 0.0_dp   ! RHS siempre 0 -> T_0^{n+1} = 0

    ! -- Filas interiores i = 2..N (Python i = 1..N-1)
    DO i = 2, N
      A(i, i-1) = -r / 2.0_dp
      A(i, i  ) =  1.0_dp + r
      A(i, i+1) = -r / 2.0_dp

      B(i, i-1) =  r / 2.0_dp
      B(i, i  ) =  1.0_dp - r
      B(i, i+1) =  r / 2.0_dp
    END DO

    ! -- Fila SZ (Python fila N): Neumann (dT/dx = 0 en x = l)
    !   Sustituyendo T_{N+1} = (4*T_N - T_{N-1})/3 en CN:
    !   Lado implicito:
    !     -r/2*T_{N-1}^{n+1} + (1+r)*T_N^{n+1} - r/2*(4T_N-T_{N-1})/3
    !     = -r/3*T_{N-1}^{n+1} + (1+r/3)*T_N^{n+1}
    A(SZ, SZ-1) = -r / 3.0_dp
    A(SZ, SZ  ) =  1.0_dp + r / 3.0_dp

    ! Lado explicito (analogo):
    B(SZ, SZ-1) =  r / 3.0_dp
    B(SZ, SZ  ) =  1.0_dp - r / 3.0_dp

  END SUBROUTINE construir_matrices


  ! ----------------------------------------------------------
  !  SUBRUTINA: imprimir_5x5
  !  Imprime las primeras 5x5 de una matriz densa.
  !  Equivalente a: print(np.round(A[:5, :5], 4))
  ! ----------------------------------------------------------
  SUBROUTINE imprimir_5x5(A)
    REAL(dp), INTENT(IN) :: A(SZ, SZ)
    INTEGER :: i, j, lim

    lim = MIN(SZ, 5)
    DO i = 1, lim
      WRITE(*, '(2X,A)', ADVANCE='NO') '['
      DO j = 1, lim
        WRITE(*, '(F7.4,A)', ADVANCE='NO') A(i, j), ' '
      END DO
      WRITE(*, '(A)') ']'
    END DO
    WRITE(*,*)

  END SUBROUTINE imprimir_5x5


  ! ----------------------------------------------------------
  !  SUBRUTINA: lu_factor_sub
  !  Factorizacion LU con pivoteo parcial. Equivalente a
  !  scipy.linalg.lu_factor(A).
  !
  !  Recibe A (SZ x SZ), devuelve:
  !    LU  : matriz con L (sin diagonal) y U solapadas
  !    piv : vector de indices de pivoteo (1..SZ, base 1)
  !
  !  Convencion identica a scipy / C++:
  !    piv(k) = indice de la fila con la que se intercambio
  !             la fila k durante la eliminacion.
  ! ----------------------------------------------------------
  SUBROUTINE lu_factor_sub(A, LU, piv)
    REAL(dp), INTENT(IN)  :: A(SZ, SZ)
    REAL(dp), INTENT(OUT) :: LU(SZ, SZ)
    INTEGER,  INTENT(OUT) :: piv(SZ)
    INTEGER  :: k, i, j, max_idx
    REAL(dp) :: max_val, tmp_row(SZ)

    ! -- Copiar A en LU
    LU = A

    DO k = 1, SZ
      ! -- Buscar pivote maximo en la columna k (pivoteo parcial)
      max_idx = k
      max_val = ABS(LU(k, k))
      DO i = k+1, SZ
        IF (ABS(LU(i, k)) > max_val) THEN
          max_val = ABS(LU(i, k))
          max_idx = i
        END IF
      END DO
      piv(k) = max_idx

      ! -- Intercambiar filas k y max_idx
      IF (max_idx /= k) THEN
        tmp_row      = LU(k,       :)
        LU(k,       :) = LU(max_idx, :)
        LU(max_idx, :) = tmp_row
      END IF

      ! -- Eliminacion gaussiana
      IF (LU(k, k) /= 0.0_dp) THEN
        DO i = k+1, SZ
          LU(i, k) = LU(i, k) / LU(k, k)
          DO j = k+1, SZ
            LU(i, j) = LU(i, j) - LU(i, k) * LU(k, j)
          END DO
        END DO
      END IF
    END DO

  END SUBROUTINE lu_factor_sub


  ! ----------------------------------------------------------
  !  SUBRUTINA: lu_solve_sub
  !  Resuelve A*x = b dado la factorizacion LU con pivoteo.
  !  Equivalente a scipy.linalg.lu_solve((lu, piv), rhs).
  !
  !  b se pasa por valor (copia interna) y x se devuelve
  !  como INTENT(OUT), igual que Vec lu_solve(...) en C++.
  ! ----------------------------------------------------------
  SUBROUTINE lu_solve_sub(LU, piv, b, x)
    REAL(dp), INTENT(IN)  :: LU(SZ, SZ)
    INTEGER,  INTENT(IN)  :: piv(SZ)
    REAL(dp), INTENT(IN)  :: b(SZ)
    REAL(dp), INTENT(OUT) :: x(SZ)
    REAL(dp) :: b2(SZ), tmp
    INTEGER  :: k, i, j

    ! -- Copiar b en b2 para no modificar el original
    b2 = b

    ! -- Aplicar permutaciones de fila al RHS
    DO k = 1, SZ
      IF (piv(k) /= k) THEN
        tmp        = b2(k)
        b2(k)      = b2(piv(k))
        b2(piv(k)) = tmp
      END IF
    END DO

    ! -- Sustitucion hacia adelante: L * y = b
    DO i = 2, SZ
      DO j = 1, i-1
        b2(i) = b2(i) - LU(i, j) * b2(j)
      END DO
    END DO

    ! -- Sustitucion hacia atras: U * x = y
    DO i = SZ, 1, -1
      DO j = i+1, SZ
        b2(i) = b2(i) - LU(i, j) * b2(j)
      END DO
      b2(i) = b2(i) / LU(i, i)
    END DO

    x = b2

  END SUBROUTINE lu_solve_sub


  ! ----------------------------------------------------------
  !  SUBRUTINA: mat_vec
  !  Multiplica matriz densa M por vector v: res = M * v
  !  Equivalente a: rhs = B @ T_actual
  ! ----------------------------------------------------------
  SUBROUTINE mat_vec(M, v, res)
    REAL(dp), INTENT(IN)  :: M(SZ, SZ), v(SZ)
    REAL(dp), INTENT(OUT) :: res(SZ)
    INTEGER :: i, j

    res = 0.0_dp
    DO i = 1, SZ
      DO j = 1, SZ
        res(i) = res(i) + M(i, j) * v(j)
      END DO
    END DO

  END SUBROUTINE mat_vec

END MODULE rutinas


! ============================================================
!  PROGRAMA PRINCIPAL
! ============================================================
PROGRAM cn_f90

  USE params
  USE rutinas
  IMPLICIT NONE

  ! ── Matrices densas A y B ──────────────────────────────────
  REAL(dp) :: A(SZ, SZ), B(SZ, SZ)
  REAL(dp) :: LU(SZ, SZ)
  INTEGER  :: piv(SZ)

  ! ── Vectores de trabajo ────────────────────────────────────
  REAL(dp) :: x_grid(SZ)
  REAL(dp) :: T_cn(SZ), T_actual(SZ), T_nuevo(SZ), rhs(SZ)

  ! ── Almacenamiento completo: T_hist(nodo, paso) ───────────
  ! Fortran: T_hist(1..SZ, 1..NTOT)  <-> C++: T_hist[0..N][0..NPT]
  REAL(dp) :: T_hist(SZ, NTOT)
  REAL(dp) :: t_vals(NTOT)
  REAL(dp) :: T_frontera_cn(NTOT)

  ! ── Parametros fisicos ─────────────────────────────────────
  REAL(dp) :: alpha, L, dx, T0, dt, t_end, r
  INTEGER  :: n_tvals

  ! ── Variables auxiliares ───────────────────────────────────
  REAL(dp) :: cpu_ini, cpu_fin, t_ejec_cn
  INTEGER  :: i, k
  CHARACTER(LEN=12) :: etiq

  ! ==========================================================
  !  CELDA 1 - Parametros del problema
  !  (equivalente al primer bloque de codigo del notebook)
  ! ==========================================================

  ! -- Parametros del problema (identicos a MELIN)
  alpha  = 0.2_dp      ! Constante de difusion termica
  L      = 1.0_dp      ! Ancho de la pared [m]
  dx     = L / REAL(N, dp)   ! Paso espacial
  T0     = 100.0_dp    ! Temperatura inicial [gC]

  ! -- Parametros temporales
  dt     = 0.05_dp     ! Paso de tiempo [s] (mismo que MELIN h=0.05)
  t_end  = 30.0_dp     ! Tiempo final [s]
  n_tvals = NPT + 1    ! np.arange(0, t_end+dt, dt) -> NPT+1 valores

  ! -- Numero de Fourier (parametro de estabilidad)
  r = alpha**2 * dt / dx**2

  ! -- Grid espacial: x = np.linspace(0, L, N+1)
  DO i = 1, SZ
    x_grid(i) = REAL(i-1, dp) * dx
  END DO

  ! -- Impresion de parametros
  WRITE(*,'(A,F3.1)')   'alpha  = ', alpha
  WRITE(*,'(A,F3.1,A)') 'l      = ', L, ' m'
  WRITE(*,'(A,I3,A)')   'N      = ', N, '  nodos'
  WRITE(*,'(A,F3.1)')   'dx     = ', dx
  WRITE(*,'(A,F4.2)')   'dt     = ', dt
  WRITE(*,'(A,I4,A)')   'NPT    = ', NPT, ' pasos de tiempo'
  WRITE(*,'(A,F8.4)')   'r      = alpha^2*dt/dx^2 = ', r
  WRITE(*,*)
  WRITE(*,'(A)') 'Crank-Nicolson es incondicionalmente estable -> cualquier r es valido'
  WRITE(*,'(A,F6.4,A,F6.4,A)') &
    '(r = ', r, ', equivalente al MELIN: r = ', &
    alpha**2 * 0.05_dp / dx**2, ')'

  ! ==========================================================
  !  CELDA 2 - construir_matrices(N, r)
  !
  !  A, B = construir_matrices(N, r)
  !  print("Matriz A (implicita) -- primeras 5x5:")
  !  print(np.round(A[:5, :5], 4))
  !  print()
  !  print("Matriz B (explicita) -- primeras 5x5:")
  !  print(np.round(B[:5, :5], 4))
  ! ==========================================================
  CALL construir_matrices(r, A, B)

  WRITE(*,*)
  WRITE(*,'(A)') 'Matriz A (implicita) -- primeras 5x5:'
  CALL imprimir_5x5(A)

  WRITE(*,'(A)') 'Matriz B (explicita) -- primeras 5x5:'
  CALL imprimir_5x5(B)

  ! ==========================================================
  !  CELDA 3 - Condicion inicial, lu_factor, bucle temporal
  !
  !  T_cn        = np.full(N+1, T0)
  !  T_cn[0]     = 0.0
  !  t_vals      = np.arange(0, t_end + dt, dt)
  !  T_hist      = np.zeros((N+1, len(t_vals)))
  !  T_hist[:,0] = T_cn.copy()
  !  lu, piv     = lu_factor(A)          <- UNA sola vez
  !  t_inicio_cn = time.perf_counter()
  !  for k in range(1, len(t_vals)):
  !      rhs        = B @ T_actual
  !      rhs[0]     = 0.0
  !      T_nuevo    = lu_solve((lu, piv), rhs)
  !      T_nuevo[0] = 0.0
  !      T_hist[:,k]= T_nuevo
  !      T_actual   = T_nuevo
  !  t_ejec_cn = time.perf_counter() - t_inicio_cn
  !  T_frontera_cn = (4*T_hist[N,:] - T_hist[N-1,:])/3.0
  ! ==========================================================

  ! -- Condicion inicial
  ! T_cn = np.full(N+1, T0)
  T_cn    = T0
  T_cn(1) = 0.0_dp   ! Dirichlet T(0,0) = 0

  ! -- Almacenamiento de la solucion
  ! t_vals = np.arange(0, t_end + dt, dt)   -> NPT+1 valores
  DO k = 1, n_tvals
    t_vals(k) = REAL(k-1, dp) * dt
  END DO

  ! T_hist = np.zeros((N+1, len(t_vals)))
  T_hist = 0.0_dp

  ! T_hist[:, 0] = T_cn.copy()
  T_hist(:, 1) = T_cn

  ! -- Factorizacion LU de A (solo una vez, fuera del bucle)
  ! lu, piv = lu_factor(A)
  CALL lu_factor_sub(A, LU, piv)

  ! -- Cronometro: t_inicio_cn = time.perf_counter()
  CALL CPU_TIME(cpu_ini)

  ! -- Bucle temporal
  ! T_actual = T_cn.copy()
  T_actual = T_cn

  DO k = 2, n_tvals

    ! rhs = B @ T_actual
    CALL mat_vec(B, T_actual, rhs)

    ! rhs[0] = 0.0   (Dirichlet: T_0 siempre 0)
    rhs(1) = 0.0_dp

    ! T_nuevo = lu_solve((lu, piv), rhs)
    CALL lu_solve_sub(LU, piv, rhs, T_nuevo)

    ! T_nuevo[0] = 0.0   (reforzar Dirichlet)
    T_nuevo(1) = 0.0_dp

    ! T_hist[:, k] = T_nuevo
    T_hist(:, k) = T_nuevo

    ! T_actual = T_nuevo
    T_actual = T_nuevo

  END DO

  ! -- t_ejec_cn = time.perf_counter() - t_inicio_cn
  CALL CPU_TIME(cpu_fin)
  t_ejec_cn = cpu_fin - cpu_ini

  ! -- Temperatura en la frontera aislada (x=l)
  ! T_frontera_cn = (4*T_hist[N, :] - T_hist[N-1, :]) / 3.0
  ! En Fortran: nodo N de Python = indice SZ en Fortran
  !             nodo N-1 de Python = indice SZ-1 en Fortran
  DO k = 1, n_tvals
    T_frontera_cn(k) = (4.0_dp*T_hist(SZ, k) - T_hist(SZ-1, k)) / 3.0_dp
  END DO

  ! -- Impresion del resultado de la celda 3
  WRITE(*,'(A,I4,A)') 'Integracion completada: ', n_tvals, ' pasos'
  WRITE(*,'(A)') '========================================='
  WRITE(*,'(A,F12.6,A)') '  Tiempo de ejecucion : ', t_ejec_cn, ' segundos'
  WRITE(*,'(A)') '  Archivo generado    : CNF.dat'
  WRITE(*,'(A)') '========================================='

  ! ==========================================================
  !  CELDA 4 - Guardar CNPY.dat + vista previa de 5 filas
  !
  !  datos_cn = np.column_stack([
  !      t_vals,
  !      T_hist[:N, :].T,    ! nodos Python 0..N-1 (incluye T0=0)
  !      T_frontera_cn
  !  ])
  !  header = "t      T1 ... T10  T_frontera"
  !  np.savetxt('CNPY.dat', datos_cn, fmt='%10.4f', ...)
  !
  !  print("Primeras 5 filas de CNPY.dat:")
  !  display(df_preview.round(4))
  !
  !  Nota de indices:
  !    Python T_hist[:N, :].T = nodos Python 0..N-1
  !                           = Fortran T_hist(1..N, k)
  !    La primera columna de datos es T_hist(1,k) = 0 (Dirichlet)
  ! ==========================================================

  OPEN(UNIT=2, FILE='CNF.dat', STATUS='REPLACE', ACTION='WRITE')

  ! -- Encabezado: t | T1 | T2 | ... | T10 | T_frontera
  WRITE(2, '(A10)', ADVANCE='NO') 't'
  DO i = 1, N
    WRITE(etiq, '(A1,I0)') 'T', i
    WRITE(2, '(A10)', ADVANCE='NO') TRIM(etiq)
  END DO
  WRITE(2, '(A12)') 'T_frontera'

  ! -- Datos: t | T_hist(1..N, k) | T_frontera_cn(k)
  !   T_hist(1..N, k) = nodos Python 0..N-1 (fmt='%10.4f')
  DO k = 1, n_tvals
    WRITE(2, '(F10.4)', ADVANCE='NO') t_vals(k)
    DO i = 1, N
      WRITE(2, '(F10.4)', ADVANCE='NO') T_hist(i, k)
    END DO
    WRITE(2, '(F12.4)') T_frontera_cn(k)
  END DO

  CLOSE(2)

  ! -- Vista previa de las primeras 5 filas
  WRITE(*,*)
  WRITE(*,'(A)') 'Primeras 5 filas de CNF.dat:'
  WRITE(*, '(A10)', ADVANCE='NO') 't'
  DO i = 1, N
    WRITE(etiq, '(A1,I0)') 'T', i
    WRITE(*, '(A10)', ADVANCE='NO') TRIM(etiq)
  END DO
  WRITE(*, '(A12)') 'T_frontera'
  WRITE(*,'(A)') REPEAT('-', 10 + 10*N + 12)

  DO k = 1, 5
    WRITE(*, '(F10.4)', ADVANCE='NO') t_vals(k)
    DO i = 1, N
      WRITE(*, '(F10.4)', ADVANCE='NO') T_hist(i, k)
    END DO
    WRITE(*, '(F12.4)') T_frontera_cn(k)
  END DO

END PROGRAM cn_f90