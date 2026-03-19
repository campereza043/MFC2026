! ============================================================
!  Solucion de la Ecuacion de Difusion
!  Metodo de Diferencias Finitas Explicito
!  Traduccion literal del notebook: EXP.ipynb  (via EXP.cpp)
!
!  CELDA 1 - Parametros del problema
!  CELDA 2 - construir_matriz_explicita(r, B) : matriz densa B
!  CELDA 3 - Condicion inicial, bucle temporal T^{n+1}=B*T^n,
!             T_frontera, tiempo de ejecucion
!  CELDA 4 - Guardar FTCSF.dat, vista previa de 5 filas
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
!  Esquema EXPLICITO:
!    T^{n+1} = B * T^n      (solo un producto matriz-vector)
!    No hay sistema lineal que resolver.
!
!    Nodos interiores:  T_i^{n+1} = r*T_{i-1} + (1-2r)*T_i + r*T_{i+1}
!    Frontera Neumann:  T_N^{n+1} = (2r/3)*T_{N-1} + (1-2r/3)*T_N
!                       (usando T_{N+1} = (4T_N - T_{N-1})/3)
!
!  Condicion de estabilidad: r = alpha^2*dt/dx^2 <= 0.5
!
!  Salida:
!    FTCSF.dat  ->  t | T1 | T2 | ... | T10 | T_frontera
!    Consola  ->  parametros, matriz 5x5, tiempo, 5 filas
!
!  Compilar:
!    gfortran -O2 -o exp_f90 EXP.f90
!  Ejecutar:
!    ./exp_f90
! ============================================================


! ============================================================
!  MODULO: params_exp
!  Centraliza precision y dimensiones.
! ============================================================
MODULE params_exp
  IMPLICIT NONE
  INTEGER,  PARAMETER :: dp   = KIND(1.0D0)   ! precision doble
  INTEGER,  PARAMETER :: N    = 10             ! divisiones espaciales
  INTEGER,  PARAMETER :: NPT  = 600            ! pasos de tiempo
  INTEGER,  PARAMETER :: SZ   = N + 1          ! tamano del sistema (0..N)
  INTEGER,  PARAMETER :: NTOT = NPT + 1        ! tiempos totales  (0..NPT)
END MODULE params_exp


! ============================================================
!  MODULO: rutinas_exp
!  Contiene todas las subrutinas en el mismo orden que el C++:
!    - construir_matriz_explicita   (CELDA 2)
!    - imprimir_5x5                 (equivalente a print numpy)
!    - mat_vec                      (equivalente a B @ T_actual)
! ============================================================
MODULE rutinas_exp
  USE params_exp
  IMPLICIT NONE

CONTAINS

  ! ----------------------------------------------------------
  !  SUBRUTINA: construir_matriz_explicita
  !  Construye B como matriz DENSA de tamano SZ x SZ.
  !  El paso temporal es:  T^{n+1} = B * T^n
  !
  !  Equivalente exacto de:
  !    B = np.zeros((size, size))
  !    B[0, 0] = 0.0
  !    for i in range(1, N):
  !        B[i,i-1]=r; B[i,i]=1-2r; B[i,i+1]=r
  !    B[N,N-1]=2r/3; B[N,N]=1-2r/3
  !
  !  Nota de indices: Fortran 1..SZ <-> Python/C++ 0..N
  !    Fila 1    (Python 0)       : Dirichlet
  !    Filas 2..N (Python 1..N-1) : esquema explicito estandar
  !    Fila SZ   (Python N)       : Neumann
  ! ----------------------------------------------------------
  SUBROUTINE construir_matriz_explicita(r, B)
    REAL(dp), INTENT(IN)  :: r
    REAL(dp), INTENT(OUT) :: B(SZ, SZ)
    INTEGER :: i

    ! -- B = np.zeros((size, size))
    B = 0.0_dp

    ! -- Fila 1 (Python fila 0): Dirichlet T(0,t) = 0
    B(1, 1) = 0.0_dp   ! T_0^{n+1} = 0 siempre

    ! -- Filas interiores i = 2..N (Python i = 1..N-1)
    DO i = 2, N
      B(i, i-1) =  r
      B(i, i  ) =  1.0_dp - 2.0_dp*r
      B(i, i+1) =  r
    END DO

    ! -- Fila SZ (Python fila N): Neumann (dT/dx = 0 en x = l)
    !   Sustituyendo T_{N+1} = (4*T_N - T_{N-1})/3 en el esquema:
    !     T_N^{n+1} = r*T_{N-1} + (1-2r)*T_N + r*(4T_N-T_{N-1})/3
    !               = (r - r/3)*T_{N-1} + (1-2r + 4r/3)*T_N
    !               = (2r/3)*T_{N-1}   + (1 - 2r/3)*T_N
    B(SZ, SZ-1) =  2.0_dp*r / 3.0_dp
    B(SZ, SZ  ) =  1.0_dp - 2.0_dp*r / 3.0_dp

  END SUBROUTINE construir_matriz_explicita


  ! ----------------------------------------------------------
  !  SUBRUTINA: imprimir_5x5
  !  Imprime las primeras 5x5 de una matriz densa.
  !  Equivalente a: print(np.round(B[:5, :5], 4))
  ! ----------------------------------------------------------
  SUBROUTINE imprimir_5x5(B)
    REAL(dp), INTENT(IN) :: B(SZ, SZ)
    INTEGER :: i, j, lim

    lim = MIN(SZ, 5)
    DO i = 1, lim
      WRITE(*, '(2X,A)', ADVANCE='NO') '['
      DO j = 1, lim
        WRITE(*, '(F7.4,A)', ADVANCE='NO') B(i, j), ' '
      END DO
      WRITE(*, '(A)') ']'
    END DO
    WRITE(*,*)

  END SUBROUTINE imprimir_5x5


  ! ----------------------------------------------------------
  !  SUBRUTINA: mat_vec
  !  Multiplica matriz densa M por vector v: res = M * v
  !  Equivalente a: T_nuevo = B @ T_actual
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

END MODULE rutinas_exp


! ============================================================
!  PROGRAMA PRINCIPAL
! ============================================================
PROGRAM exp_f90

  USE params_exp
  USE rutinas_exp
  IMPLICIT NONE

  ! ── Matriz densa B ─────────────────────────────────────────
  REAL(dp) :: B(SZ, SZ)

  ! ── Vectores de trabajo ────────────────────────────────────
  REAL(dp) :: x_grid(SZ)
  REAL(dp) :: T_exp(SZ), T_actual(SZ), T_nuevo(SZ)

  ! ── Almacenamiento completo: T_hist(nodo, paso) ───────────
  ! Fortran: T_hist(1..SZ, 1..NTOT)  <-> C++: T_hist[0..N][0..NPT]
  REAL(dp) :: T_hist(SZ, NTOT)
  REAL(dp) :: t_vals(NTOT)
  REAL(dp) :: T_frontera(NTOT)

  ! ── Parametros fisicos ─────────────────────────────────────
  REAL(dp) :: alpha, L, dx, T0, dt, t_end, r

  ! ── Variables auxiliares ───────────────────────────────────
  REAL(dp) :: cpu_ini, cpu_fin, t_ejec
  INTEGER  :: i, k, n_tvals
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
  IF (r <= 0.5_dp) THEN
    WRITE(*,'(A,F6.4,A)') 'Criterio de estabilidad: r = ', r, ' <= 0.5 -> ESTABLE'
  ELSE
    WRITE(*,'(A,F6.4,A)') 'Criterio de estabilidad: r = ', r, ' >  0.5 -> INESTABLE'
  END IF

  ! ==========================================================
  !  CELDA 2 - construir_matriz_explicita(r, B)
  !
  !  B = construir_matriz_explicita(N, r)
  !  print("Matriz B (explicita) -- primeras 5x5:")
  !  print(np.round(B[:5, :5], 4))
  ! ==========================================================
  CALL construir_matriz_explicita(r, B)

  WRITE(*,*)
  WRITE(*,'(A)') 'Matriz B (explicita) -- primeras 5x5:'
  CALL imprimir_5x5(B)

  ! ==========================================================
  !  CELDA 3 - Condicion inicial, bucle temporal
  !
  !  T_exp       = np.full(N+1, T0)
  !  T_exp[0]    = 0.0
  !  t_vals      = np.arange(0, t_end + dt, dt)
  !  T_hist      = np.zeros((N+1, len(t_vals)))
  !  T_hist[:,0] = T_exp.copy()
  !  t_inicio    = time.perf_counter()
  !  for k in range(1, len(t_vals)):
  !      T_nuevo    = B @ T_actual    <- avance explicito
  !      T_nuevo[0] = 0.0
  !      T_hist[:,k]= T_nuevo
  !      T_actual   = T_nuevo
  !  t_ejec = time.perf_counter() - t_inicio
  !  T_frontera = (4*T_hist[N,:] - T_hist[N-1,:])/3.0
  ! ==========================================================

  ! -- Condicion inicial
  ! T_exp = np.full(N+1, T0)
  T_exp    = T0
  T_exp(1) = 0.0_dp   ! Dirichlet T(0,0) = 0

  ! -- Almacenamiento de la solucion
  ! t_vals = np.arange(0, t_end + dt, dt)   -> NPT+1 valores
  DO k = 1, n_tvals
    t_vals(k) = REAL(k-1, dp) * dt
  END DO

  ! T_hist = np.zeros((N+1, len(t_vals)))
  T_hist = 0.0_dp

  ! T_hist[:, 0] = T_exp.copy()
  T_hist(:, 1) = T_exp

  ! -- Cronometro: t_inicio = time.perf_counter()
  CALL CPU_TIME(cpu_ini)

  ! -- Bucle temporal: T^{n+1} = B * T^n
  ! T_actual = T_exp.copy()
  T_actual = T_exp

  DO k = 2, n_tvals

    ! T_nuevo = B @ T_actual   (avance explicito)
    CALL mat_vec(B, T_actual, T_nuevo)

    ! T_nuevo[0] = 0.0   (reforzar Dirichlet)
    T_nuevo(1) = 0.0_dp

    ! T_hist[:, k] = T_nuevo
    T_hist(:, k) = T_nuevo

    ! T_actual = T_nuevo
    T_actual = T_nuevo

  END DO

  ! -- t_ejec = time.perf_counter() - t_inicio
  CALL CPU_TIME(cpu_fin)
  t_ejec = cpu_fin - cpu_ini

  ! -- Temperatura en la frontera aislada (x=l)
  ! T_frontera = (4*T_hist[N, :] - T_hist[N-1, :]) / 3.0
  ! En Fortran: nodo N de Python = indice SZ en Fortran
  !             nodo N-1 de Python = indice SZ-1 en Fortran
  DO k = 1, n_tvals
    T_frontera(k) = (4.0_dp*T_hist(SZ, k) - T_hist(SZ-1, k)) / 3.0_dp
  END DO

  ! -- Impresion del resultado de la celda 3
  WRITE(*,'(A,I4,A)') 'Integracion completada: ', n_tvals, ' pasos'
  WRITE(*,'(A)') '========================================='
  WRITE(*,'(A,F12.6,A)') '  Tiempo de ejecucion : ', t_ejec, ' segundos'
  WRITE(*,'(A)') '  Archivo generado    : FTCSF.dat'
  WRITE(*,'(A)') '========================================='

  ! ==========================================================
  !  CELDA 4 - Guardar FTCSF.dat + vista previa de 5 filas
  !
  !  datos_exp = np.column_stack([
  !      t_vals,
  !      T_hist[:N, :].T,   ! nodos Python 0..N-1 (T0=0 incluido)
  !      T_frontera
  !  ])
  !  header = "t      T1 ... T10  T_frontera"
  !  np.savetxt('FTCSF.dat', datos_exp, fmt='%10.4f', ...)
  !
  !  Nota de indices:
  !    Python T_hist[:N, :].T = nodos Python 0..N-1
  !                           = Fortran T_hist(1..N, k)
  !    Primera columna de datos = T_hist(1,k) = 0 (Dirichlet)
  ! ==========================================================

  OPEN(UNIT=2, FILE='FTCSF.dat', STATUS='REPLACE', ACTION='WRITE')

  ! -- Encabezado: t | T1 | T2 | ... | T10 | T_frontera
  WRITE(2, '(A10)', ADVANCE='NO') 't'
  DO i = 1, N
    WRITE(etiq, '(A1,I0)') 'T', i
    WRITE(2, '(A10)', ADVANCE='NO') TRIM(etiq)
  END DO
  WRITE(2, '(A12)') 'T_frontera'

  ! -- Datos: t | T_hist(1..N, k) | T_frontera(k)
  !   T_hist(1..N, k) = nodos Python 0..N-1 (fmt='%10.4f')
  DO k = 1, n_tvals
    WRITE(2, '(F10.4)', ADVANCE='NO') t_vals(k)
    DO i = 1, N
      WRITE(2, '(F10.4)', ADVANCE='NO') T_hist(i, k)
    END DO
    WRITE(2, '(F12.4)') T_frontera(k)
  END DO

  CLOSE(2)

  ! -- Vista previa de las primeras 5 filas
  WRITE(*,*)
  WRITE(*,'(A)') 'Primeras 5 filas de FTCSF.dat:'
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
    WRITE(*, '(F12.4)') T_frontera(k)
  END DO

END PROGRAM exp_f90