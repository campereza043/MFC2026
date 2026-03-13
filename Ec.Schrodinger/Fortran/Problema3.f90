program schrodinger_cn_2
    implicit none

    integer, parameter    :: m       = 20
    real(8), parameter    :: k       = 0.0014d0
    real(8), parameter    :: pi      = 3.141592653589793d0
    real(8), parameter    :: a       = 0.0d0
    real(8), parameter    :: b       = pi
    real(8), parameter    :: T_final = 1.0d0
    integer, parameter    :: N_steps = nint(T_final / k)
    complex(8), parameter :: eye     = (0.0d0, 1.0d0)

    integer     :: ii, nn, n_int
    real(8)     :: h, lam, t_n, t_np1, t_final_efectivo
    real(8)     :: start_time, end_time, suma_error_sq, error_global

    complex(8), allocatable :: U(:), dl(:), d(:), du(:), rhs(:), sol(:)
    complex(8) :: bc_L_n, bc_R_n, bc_L_np1, bc_R_np1, u_analitica

    call cpu_time(start_time)

    h     = (b - a) / dble(m)
    lam   = k / (h**2)
    n_int = m - 1

    print *, "--- Simulación: i*u_t = -u_xx | u = exp(x + it) ---"
    print '(A,F10.8,A,F8.6,A,I6)', " Delta x: ", h, " | Delta t: ", k, " | Pasos: ", N_steps

    allocate(U(0:m), dl(n_int), d(n_int), du(n_int), rhs(n_int), sol(n_int))

    ! Condición inicial: u(x,0) = e^x
    do ii = 0, m
        U(ii) = exp(a + ii*h)
    end do

    ! Diagonales de A (igual que C++: diag_A = 2i - lam, fuera = +lam/2)
    dl = lam / 2.0d0          ! subdiagonal  de A: +lam/2
    d  = 2.0d0*eye - lam      ! diagonal     de A:  2i - lam
    du = lam / 2.0d0          ! superdiagonal de A: +lam/2

    open(unit=20, file='raro_1000.dat', status='replace')
    write(20, '(A)') "# t            x            u_real       u_imag       u_abs"

    t_final_efectivo = 0.0d0

    do nn = 0, N_steps - 1
        t_n   = nn * k
        t_np1 = (nn + 1) * k
        t_final_efectivo = t_np1

        bc_L_n   = exp(a + eye * t_n)
        bc_R_n   = exp(b + eye * t_n)
        bc_L_np1 = exp(a + eye * t_np1)
        bc_R_np1 = exp(b + eye * t_np1)

        ! rhs = B * U_interior  (B: diag = 2i+lam, fuera = -lam/2)
        do ii = 1, n_int
            rhs(ii) = (2.0d0*eye + lam) * U(ii)
            if (ii > 1)     rhs(ii) = rhs(ii) - (lam/2.0d0) * U(ii-1)
            if (ii < n_int) rhs(ii) = rhs(ii) - (lam/2.0d0) * U(ii+1)
        end do

        ! Contribución de frontera: C(0) = -(lam/2)*(bc_L_n + bc_L_np1)
        ! pero A tiene +lam/2 fuera, así que el término que pasa al rhs es:
        ! rhs(1)     -= (+lam/2)*bc_L_np1  y  rhs(1) ya tiene el de B
        ! Replicando exactamente el C++:
        rhs(1)     = rhs(1)     - (lam/2.0d0) * (bc_L_n + bc_L_np1)
        rhs(n_int) = rhs(n_int) - (lam/2.0d0) * (bc_R_n + bc_R_np1)

        call thomas_complex(dl, d, du, rhs, sol, n_int)

        U(0)     = bc_L_np1
        U(m)     = bc_R_np1
        U(1:m-1) = sol

        do ii = 0, m
            write(20, '(5ES16.6)') &
                t_np1, a + ii*h, real(U(ii)), aimag(U(ii)), abs(U(ii))
        end do
        write(20, *)
    end do

    close(20)
    call cpu_time(end_time)

    suma_error_sq = 0.0d0
    do ii = 0, m
        u_analitica   = exp((a + ii*h) + eye * t_final_efectivo)
        suma_error_sq = suma_error_sq + abs(U(ii) - u_analitica)**2
    end do
    error_global = sqrt(suma_error_sq)

    print *, "------------------------------------"
    print *, "Finalizado en Debian."
    print '(A, F10.6, A)', " Tiempo de ejecución: ", (end_time - start_time), " s"
    print '(A, E13.6)',    " Error Global (Norma L2): ", error_global
    print *, "------------------------------------"

    deallocate(U, dl, d, du, rhs, sol)

contains

    subroutine thomas_complex(dl_v, d_v, du_v, b_vec, x_vec, n_size)
        integer,    intent(in)  :: n_size
        complex(8), dimension(n_size), intent(in)  :: dl_v, d_v, du_v, b_vec
        complex(8), dimension(n_size), intent(out) :: x_vec
        complex(8), dimension(n_size) :: cp, dp
        integer :: i

        cp(1) = du_v(1) / d_v(1)
        dp(1) = b_vec(1) / d_v(1)

        do i = 2, n_size
            cp(i) = du_v(i) / (d_v(i) - dl_v(i)*cp(i-1))
            dp(i) = (b_vec(i) - dl_v(i)*dp(i-1)) / (d_v(i) - dl_v(i)*cp(i-1))
        end do

        x_vec(n_size) = dp(n_size)
        do i = n_size-1, 1, -1
            x_vec(i) = dp(i) - cp(i)*x_vec(i+1)
        end do
    end subroutine thomas_complex

end program schrodinger_cn_2