program schrodinger_cn
    implicit none

    ! Parámetros del problema
    integer, parameter :: m = 1000
    real(8), parameter :: k = 0.0013d0
    real(8), parameter :: a = -2.0d0
    real(8), parameter :: b = 2.0d0
    integer, parameter :: N_steps = 769 
    complex(8), parameter :: eye = (0.0d0, 1.0d0)

    ! Variables
    integer :: ii, nn, n_int
    real(8) :: h, lam, t_actual, t_np1, start_time, end_time
    real(8), dimension(0:m) :: x_full
    complex(8), dimension(0:N_steps, 0:m) :: U_full
    complex(8), dimension(m-1) :: rhs, sol, dl, d, du
    real(8) :: suma_cuadrados, error_l2, err_local

    call cpu_time(start_time)

    h = (b - a) / dble(m)
    lam = k / (h**2)
    n_int = m - 1

    print *, "h =", h
    print *, "lambda =", lam

    ! Condición Inicial
    do ii = 0, m
        x_full(ii) = a + ii * h
        U_full(0, ii) = phi(x_full(ii))
    end do

    ! Evolución Temporal
    do nn = 0, N_steps - 1
        t_actual = nn * k
        t_np1 = (nn + 1) * k

        do ii = 1, n_int
            rhs(ii) = (2.0d0*eye - 2.0d0*lam) * U_full(nn, ii)
            if (ii > 1) rhs(ii) = rhs(ii) + lam * U_full(nn, ii-1)
            if (ii < n_int) rhs(ii) = rhs(ii) + lam * U_full(nn, ii+1)
        end do

        rhs(1) = rhs(1) + lam * (bc_left(t_actual) + bc_left(t_np1))
        rhs(n_int) = rhs(n_int) + lam * (bc_right(t_actual) + bc_right(t_np1))

        dl = -lam
        d  = 2.0d0*eye + 2.0d0*lam
        du = -lam
        
        call thomas_complex(dl, d, du, rhs, sol, n_int)

        U_full(nn+1, 0) = bc_left(t_np1)
        U_full(nn+1, m) = bc_right(t_np1)
        U_full(nn+1, 1:m-1) = sol
    end do

    call cpu_time(end_time)
    
    ! Cálculo de error L2 global en el tiempo final
    suma_cuadrados = 0.0d0
    do ii = 0, m
        err_local = abs(U_full(N_steps, ii) - u_exacta(x_full(ii), N_steps*k))
        suma_cuadrados = suma_cuadrados + err_local**2
    end do
    
    error_l2 = sqrt(suma_cuadrados / dble(m + 1))

    print '(A, E15.7)', " Error L2 global final =", error_l2
    print '(A, F10.7, A)', " Tiempo de ejecucion =", (end_time - start_time), " segundos"

    ! EXPORTAR A .DAT
    print *, "Guardando resultados en resultados.dat..."
    call export_dat(U_full, x_full, N_steps, m, k)
    print *, "¡Archivo .dat generado con éxito!"

contains

    complex(8) function phi(x)
        real(8), intent(in) :: x
        phi = cmplx(sin(x), 0.0d0, 8)
    end function

    complex(8) function bc_left(t)
        real(8), intent(in) :: t
        bc_left = exp(eye * t) * sin(-2.0d0)
    end function

    complex(8) function bc_right(t)
        real(8), intent(in) :: t
        bc_right = exp(eye * t) * sin(2.0d0)
    end function

    complex(8) function u_exacta(x, t)
        real(8), intent(in) :: x, t
        u_exacta = exp(eye * t) * sin(x)
    end function

    subroutine thomas_complex(dl, d, du, b_vec, x_vec, n_size)
        integer, intent(in) :: n_size
        complex(8), dimension(n_size), intent(in) :: dl, d, du, b_vec
        complex(8), dimension(n_size), intent(out) :: x_vec
        complex(8), dimension(n_size) :: cp, dp
        integer :: i
        cp(1) = du(1) / d(1)
        dp(1) = b_vec(1) / d(1)
        do i = 2, n_size
            cp(i) = du(i) / (d(i) - dl(i)*cp(i-1))
            dp(i) = (b_vec(i) - dl(i)*dp(i-1)) / (d(i) - dl(i)*cp(i-1))
        end do
        x_vec(n_size) = dp(n_size)
        do i = n_size-1, 1, -1
            x_vec(i) = dp(i) - cp(i)*x_vec(i+1)
        end do
    end subroutine thomas_complex

    subroutine export_dat(U, x, steps, space_m, dt)
        complex(8), dimension(0:steps, 0:space_m), intent(in) :: U
        real(8), dimension(0:space_m), intent(in) :: x
        integer, intent(in) :: steps, space_m
        real(8), intent(in) :: dt
        integer :: t_idx, x_idx
        
        open(unit=20, file='resultados.dat', status='replace')
        write(20, '(A)') "# t            x            u_real       u_imag       u_abs"
        
        do t_idx = 0, steps
            do x_idx = 0, space_m
                write(20, '(5F13.6)') &
                    t_idx*dt, x(x_idx), real(U(t_idx, x_idx)), &
                    aimag(U(t_idx, x_idx)), abs(U(t_idx, x_idx))
            end do
            write(20, *) 
        end do
        close(20)
    end subroutine export_dat

end program schrodinger_cn