! Resuelve la ecuación estacionaria de Schrödinger para el potencial V(x) = x⁴
! utilizando el método de diferencias finitas de segundo orden.
PROGRAM Osciladorx4
    implicit none
    integer :: i, j, k, N, Nf
    real(8), allocatable :: A(:,:), d(:), e(:), tem(:)
    real(8) :: h, xi, h_2, Rmin, Rmax, aux_val
    real(8), allocatable :: aux_vec(:)

    ! Configuración del dominio
    Rmin = -10.d0
    Rmax = 10.d0
    Nf = 1600 ! Puntos máximos en la malla

    write(6,*) 'Calculando los 6 primeros autovalores para V(x) = x⁴:'
    write(6,*) '   N      E0        E1        E2        E3        E4        E5'
    write(6,*) '------------------------------------------------------------'

    N = 50
    do while(N .le. Nf)
        h = (Rmax - Rmin) / dble(N)
        h_2 = 1.d0 / h**2
        
        ! Asignación de memoria
        allocate(A(N,N), d(N), e(N), tem(N), aux_vec(N))
        
        A = 0.d0
        d = 0.d0
        e = 0.d0
        
        ! Construcción de la matriz Hamiltoniana (Diferencias Finitas)
        ! Representamos H = -1/2 * d2/dx2 + |x|
        ! Multiplicamos por 2 para diagonalizar: 2H = -d2/dx2 + 2|x|
        do i = 1, N
            xi = Rmin + i * h
            d(i) = 2.d0 * h_2 + 2.d0 * xi**4 ! Diagonal: 2/h^2 + 2*V(x)
            e(i) = -h_2                        ! Extradiagonal: -1/h^2
            A(i,i) = 1.d0                      ! Matriz de autovectores inicial
        end do

        ! Subrutina de diagonalización (debe estar disponible en tu librería)
        call diagotri(d, e, N, A, .true.)

        ! Ordenamiento de autovalores y autovectores (Bubble Sort)
        do i = 1, N - 1
            do j = i + 1, N
                if (d(j) .lt. d(i)) then
                    ! Intercambiar autovalores
                    aux_val = d(i); d(i) = d(j); d(j) = aux_val
                    ! Intercambiar autovectores (columnas de A)
                    aux_vec(:) = A(:,i)
                    A(:,i) = A(:,j)
                    A(:,j) = aux_vec(:)
                end if
            end do
        end do
        
        ! Los autovalores reales son 0.5 * d (porque multiplicamos la ec. por 2)
        write(6,"(i4, 6(2x, F9.6))") N, 0.5d0 * d(1:6)

        ! Guardado de datos en la última iteración (N = Nf)
        if (N .eq. Nf) then
            open(unit=1, file="data_x4_fortran")
            write(1,"('#', i4, 6x, 6(F10.6, 1x))") N, 0.5d0 * d(1:6)
            
            ! Cálculo de la norma para cada uno de los primeros 6 autovectores
            ! Norma = Integral(|psi|^2 dx)
            do j = 1, 6
                tem(j) = 0.d0
                do i = 1, N
                    tem(j) = tem(j) + A(i,j)**2
                end do
                tem(j) = tem(j) * h
            end do

            ! Guardado de xi y las densidades de probabilidad normalizadas
            do i = 1, N
                xi = Rmin + i * h
                write(1,"(F10.5, 6(1x, F10.5))") xi, &
                    A(i,1)**2/tem(1), A(i,2)**2/tem(2), A(i,3)**2/tem(3), &
                    A(i,4)**2/tem(4), A(i,5)**2/tem(5), A(i,6)**2/tem(6)
            end do
            close(1)
        end if

        deallocate(A, d, e, tem, aux_vec)
        N = N * 2
        
    end do

    write(6,*) '------------------------------------------------------------'
    write(6,*) 'Proceso finalizado. Autovectores en: data_x4_fortran'
    
END PROGRAM Osciladorx4