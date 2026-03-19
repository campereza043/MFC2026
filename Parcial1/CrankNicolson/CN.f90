! ============================================================
! Solucion de la Ecuacion de Difusion - Metodo de Crank-Nicolson
! Traduccion a Fortran 90 (Formato Libre)
! ============================================================

PROGRAM crank_nicolson
    IMPLICIT NONE
    
    ! --- Definicion de Precision (Doble Precision) ---
    INTEGER, PARAMETER :: dp = KIND(1.0d0)
    
    ! --- Parametros del Problema ---
    INTEGER, PARAMETER :: N = 10
    INTEGER, PARAMETER :: NPT = 600
    INTEGER, PARAMETER :: nmax = N + 1
    INTEGER, PARAMETER :: npmax = NPT + 1
    
    REAL(dp), PARAMETER :: alpha = 0.2_dp
    REAL(dp), PARAMETER :: ancho = 1.0_dp
    REAL(dp), PARAMETER :: t0    = 100.0_dp
    REAL(dp), PARAMETER :: dt    = 0.05_dp
    REAL(dp), PARAMETER :: t_end = 30.0_dp
    REAL(dp) :: dx, r, pi, cpu_ini, cpu_fin, t_ejec_cn
    
    ! --- Arreglos del Sistema ---
    REAL(dp) :: aa(nmax), ba(nmax), ca(nmax) ! Diagonales de A (implicita)
    REAL(dp) :: ab(nmax), bb(nmax), cb(nmax) ! Diagonales de B (explicita)
    REAL(dp) :: x(nmax), t_actual(nmax), t_nuevo(nmax), rhs(nmax)
    REAL(dp) :: t_hist(nmax, npmax), t_vals(npmax), t_front(npmax)
    
    ! --- Metricas ---
    REAL(dp) :: mae_nodo(nmax), mae_global, rmse_global, max_error
    REAL(dp) :: t_an, err, diff, diffmin
    INTEGER  :: i, k, paso, idx, total, itcomp
    REAL(dp) :: t_comp_list(3) = (/ 0.05_dp, 5.0_dp, 30.0_dp /)

    ! --- Inicializacion de constantes ---
    dx = ancho / REAL(N, dp)
    pi = ACOS(-1.0_dp)
    r  = (alpha**2 * dt) / (dx**2)
    
    WRITE(*,'(A,F6.1)')   'alpha  = ', alpha
    WRITE(*,'(A,F6.1,A)') 'l      = ', ancho, ' m'
    WRITE(*,'(A,I3,A)')   'N      = ', N,     ' nodos'
    WRITE(*,'(A,F6.4)')   'dx     = ', dx
    WRITE(*,'(A,F6.4)')   'dt     = ', dt
    WRITE(*,'(A,F8.4)')   'r      = ', r
    WRITE(*,*) 'Crank-Nicolson es incondicionalmente estable.'
    WRITE(*,*)

    ! --- Grid espacial y temporal ---
    DO i = 1, nmax
        x(i) = REAL(i-1, dp) * dx
    END DO
    DO k = 1, npmax
        t_vals(k) = REAL(k-1, dp) * dt
    END DO

    ! --- Construccion de Matrices A y B ---
    aa = 0.0_dp; ba = 0.0_dp; ca = 0.0_dp
    ab = 0.0_dp; bb = 0.0_dp; cb = 0.0_dp

    ! Fila 1: Dirichlet T(0,t) = 0
    ba(1) = 1.0_dp
    bb(1) = 0.0_dp

    ! Filas interiores (2 a N)
    DO i = 2, N
        aa(i) = -r/2.0_dp; ba(i) = 1.0_dp + r; ca(i) = -r/2.0_dp
        ab(i) =  r/2.0_dp; bb(i) = 1.0_dp - r; cb(i) =  r/2.0_dp
    END DO

    ! Fila NMAX: Neumann (dT/dx = 0 en x = L)
    aa(nmax) = -r/3.0_dp
    ba(nmax) =  1.0_dp + r/3.0_dp
    ab(nmax) =  r/3.0_dp
    bb(nmax) =  1.0_dp - r/3.0_dp

    ! Imprimir Matrices 5x5
    WRITE(*,'(A)') 'Matriz A (implicita) - primeras 5x5:'
    CALL impr5x5(aa, ba, ca, nmax)
    WRITE(*,'(A)') 'Matriz B (explicita) - primeras 5x5:'
    CALL impr5x5(ab, bb, cb, nmax)

    ! --- Condicion Inicial ---
    t_actual = t0
    t_actual(1) = 0.0_dp
    t_hist(:, 1) = t_actual

    ! --- Bucle Temporal ---
    CALL CPU_TIME(cpu_ini)
    DO paso = 2, npmax
        ! RHS = B * T_actual
        CALL mult_tridiag(ab, bb, cb, t_actual, rhs, nmax)
        rhs(1) = 0.0_dp ! Condicion Dirichlet en RHS
        
        ! Resolver A * T_nuevo = RHS
        CALL thomas(aa, ba, ca, rhs, t_nuevo, nmax)
        
        t_nuevo(1) = 0.0_dp ! Refuerzo Dirichlet
        t_hist(:, paso) = t_nuevo
        t_actual = t_nuevo
    END DO
    CALL CPU_TIME(cpu_fin)
    t_ejec_cn = cpu_fin - cpu_ini

    ! Temperatura frontera aislada (O(Dx^2))
    DO k = 1, npmax
        t_front(k) = (4.0_dp * t_hist(nmax, k) - t_hist(nmax-1, k)) / 3.0_dp
    END DO

    ! --- Guardado de datos ---
    OPEN(UNIT=10, FILE='CN.dat', STATUS='REPLACE')
    WRITE(10,'(A10,10A10,A12)') 't','T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T_frontera'
    DO k = 1, npmax
        WRITE(10,'(11F10.4, F12.4)') t_vals(k), (t_hist(i,k), i=2, nmax), t_front(k)
    END DO
    CLOSE(10)

    ! --- Analisis de Error y Tablas ---
    DO itcomp = 1, 3
        ! Buscar indice mas cercano a t_comp
        idx = 1
        diffmin = ABS(t_vals(1) - t_comp_list(itcomp))
        DO k = 2, npmax
            diff = ABS(t_vals(k) - t_comp_list(itcomp))
            IF (diff < diffmin) THEN
                diffmin = diff
                idx = k
            END IF
        END DO

        WRITE(*,'(/,A,F8.4,A)') ' === t = ', t_vals(idx), ' s ==='
        WRITE(*,'(A10,A14,A16,A14)') 'x [m]', 'CN [gC]', 'Analitica', '|Error|'
        DO i = 1, nmax
            t_an = t_analitica_func(x(i), t_vals(idx), alpha, ancho, pi)
            err  = ABS(t_hist(i, idx) - t_an)
            WRITE(*,'(F10.4, F14.4, F16.4, F14.6)') x(i), t_hist(i, idx), t_an, err
        END DO
    END DO

    ! --- Metricas Globales ---
    mae_nodo = 0.0_dp; mae_global = 0.0_dp; rmse_global = 0.0_dp; max_error = 0.0_dp
    total = 0
    DO k = 1, npmax
        DO i = 1, nmax
            t_an = t_analitica_func(x(i), t_vals(k), alpha, ancho, pi)
            err = ABS(t_hist(i, k) - t_an)
            mae_nodo(i) = mae_nodo(i) + err
            mae_global  = mae_global + err
            rmse_global = rmse_global + err**2
            IF (err > max_error) max_error = err
            total = total + 1
        END DO
    END DO
    
    mae_global  = mae_global / REAL(total, dp)
    rmse_global = SQRT(rmse_global / REAL(total, dp))
    
    WRITE(*,'(/,A,F12.6)') 'MAE Global:  ', mae_global
    WRITE(*,'(A,F12.6)')   'RMSE Global: ', rmse_global
    WRITE(*,'(A,F12.6)')   'Error Max:   ', max_error

CONTAINS

    SUBROUTINE impr5x5(a, b, c, nsize)
        REAL(dp), INTENT(IN) :: a(:), b(:), c(:)
        INTEGER, INTENT(IN)  :: nsize
        INTEGER :: row, col, lim
        REAL(dp) :: val
        lim = MIN(nsize, 5)
        DO row = 1, lim
            WRITE(*,'(A)', ADVANCE='NO') ' ['
            DO col = 1, lim
                val = 0.0_dp
                IF (col == row - 1) val = a(row)
                IF (col == row)     val = b(row)
                IF (col == row + 1) val = c(row)
                WRITE(*,'(F8.4)', ADVANCE='NO') val
            END DO
            WRITE(*,'(A)') '  ]'
        END DO
        WRITE(*,*)
    END SUBROUTINE

    SUBROUTINE mult_tridiag(a, b, c, v, res, nsize)
        REAL(dp), INTENT(IN)  :: a(:), b(:), c(:), v(:)
        REAL(dp), INTENT(OUT) :: res(:)
        INTEGER, INTENT(IN)   :: nsize
        INTEGER :: j
        DO j = 1, nsize
            res(j) = b(j) * v(j)
            IF (j > 1)     res(j) = res(j) + a(j) * v(j-1)
            IF (j < nsize) res(j) = res(j) + c(j) * v(j+1)
        END DO
    END SUBROUTINE

    SUBROUTINE thomas(a, b, c, d, sol, nsize)
        REAL(dp), INTENT(IN)  :: a(:), b(:), c(:), d(:)
        REAL(dp), INTENT(OUT) :: sol(:)
        INTEGER, INTENT(IN)   :: nsize
        REAL(dp) :: c_p(nsize), d_p(nsize), m
        INTEGER :: j
        
        ! Forward elimination
        c_p(1) = c(1) / b(1)
        d_p(1) = d(1) / b(1)
        DO j = 2, nsize
            m = b(j) - a(j) * c_p(j-1)
            c_p(j) = c(j) / m
            d_p(j) = (d(j) - a(j) * d_p(j-1)) / m
        END DO
        
        ! Backward substitution
        sol(nsize) = d_p(nsize)
        DO j = nsize-1, 1, -1
            sol(j) = d_p(j) - c_p(j) * sol(j+1)
        END DO
    END SUBROUTINE

    REAL(dp) FUNCTION t_analitica_func(px, pt, palpha, pl, ppi)
        REAL(dp), INTENT(IN) :: px, pt, palpha, pl, ppi
        REAL(dp) :: suma, kn
        INTEGER :: n_idx
        suma = 0.0_dp
        DO n_idx = 0, 49
            kn = REAL(2*n_idx + 1, dp) * ppi / (2.0_dp * pl)
            suma = suma + (1.0_dp / REAL(2*n_idx + 1, dp)) * &
                   SIN(kn * px) * EXP(-palpha**2 * kn**2 * pt)
        END DO
        t_analitica_func = (400.0_dp / ppi) * suma
    END FUNCTION

END PROGRAM crank_nicolson