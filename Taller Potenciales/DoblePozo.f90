! Resuelve la ecuación estacionaria de Schrödinger para el potencial de doble pozo:
! V(x) = x^4 - A*x^2
! utilizando el método de diferencias finitas de segundo orden.
PROGRAM DoblePozo
    implicit none
    integer :: i, j, k, N, Nf
    real(8), allocatable :: A_mat(:,:), d(:), e(:), tem(:)
    real(8) :: h, xi, h_2, Rmin, Rmax, aux_val
    real(8) :: param_A
    real(8), allocatable :: aux_vec(:)

    ! --- Parámetros del Potencial ---
    param_A = 5.d0  ! Puedes variar este valor para profundizar los pozos
    ! --------------------------------

    ! Configuración del dominio
    ! Nota: Para x^4 - Ax^2, un Rmin/Rmax de -5 a 5 suele ser suficiente.
    Rmin = -5.d0
    Rmax = 5.d0
    Nf = 1600 ! Puntos máximos en la malla

    write(6,*) 'Calculando autovalores para Doble Pozo V(x) = x^4 - ', param_A, 'x^2:'
    write(6,*) '   N      E0        E1        E2        E3        E4        E5'
    write(6,*) '------------------------------------------------------------'

    N = 50
    do while(N .le. Nf)
        h = (Rmax - Rmin) / dble(N)
        h_2 = 1.d0 / h**2
        
        ! Asignación de memoria
        allocate(A_mat(N,N), d(N), e(N), tem(N), aux_vec(N))
        
        A_mat = 0.d0; d = 0.d0; e = 0.d0
        
        ! Construcción de la matriz Hamiltoniana
        ! 2H = -d2/dx2 + 2*(x^4 - A*x^2)
        do i = 1, N
            xi = Rmin + i * h
            ! DIAGONAL: 2/h^2 + 2*V(xi)
            d(i) = 2.d0 * h_2 + 2.d0 * (xi**4 - param_A * xi**2)
            e(i) = -h_2
            A_mat(i,i) = 1.d0
        end do

        ! Subrutina de diagonalización
        call diagotri(d, e, N, A_mat, .true.)

        ! Ordenamiento (Bubble Sort)
        do i = 1, N - 1
            do j = i + 1, N
                if (d(j) .lt. d(i)) then
                    aux_val = d(i); d(i) = d(j); d(j) = aux_val
                    aux_vec(:) = A_mat(:,i)
                    A_mat(:,i) = A_mat(:,j)
                    A_mat(:,j) = aux_vec(:)
                end if
            end do
        end do
        
        ! Resultados (E = 0.5 * autovalor de la matriz)
        write(6,"(i4, 6(2x, F9.6))") N, 0.5d0 * d(1:6)

        ! Guardado de datos en la última iteración
        if (N .eq. Nf) then
            open(unit=1, file="data_doble_pozo_fortran")
            write(1,"('#', i4, 6x, 6(F10.6, 1x))") N, 0.5d0 * d(1:6)
            
            ! Normalización
            do j = 1, 6
                tem(j) = 0.d0
                do i = 1, N
                    tem(j) = tem(j) + A_mat(i,j)**2
                end do
                tem(j) = tem(j) * h
            end do

            ! Guardado de densidades
            do i = 1, N
                xi = Rmin + i * h
                write(1,"(F10.5, 6(1x, F10.5))") xi, &
                    A_mat(i,1)**2/tem(1), A_mat(i,2)**2/tem(2), A_mat(i,3)**2/tem(3), &
                    A_mat(i,4)**2/tem(4), A_mat(i,5)**2/tem(5), A_mat(i,6)**2/tem(6)
            end do
            close(1)
        end if

        deallocate(A_mat, d, e, tem, aux_vec)
        N = N * 2
    end do

    write(6,*) '------------------------------------------------------------'
    write(6,*) 'Proceso finalizado. Autovectores en: data_doble_pozo_fortran'
    
END PROGRAM DoblePozo