program metodos_modelo
implicit none

! parametros
integer, parameter :: N = 50
real(8) :: rb, h, p0
real(8) :: t(0:N), p_exact(0:N)
real(8) :: p_euler(0:N), p_taylor(0:N), p_trap(0:N)
real(8) :: start, finish
real(8) :: time_euler, time_taylor, time_trap
real(8) :: error_euler, error_taylor, error_trap
integer :: i

rb = 0.002d0
h = 1.0d0
p0 = 0.01d0

! =========================
! inicializacion tiempo y solucion exacta
! =========================

do i=0,N
    t(i) = i*h
    p_exact(i) = 1.0d0 - 0.99d0*exp(-rb*t(i))
end do

! condiciones iniciales
p_euler(0) = p0
p_taylor(0) = p0
p_trap(0) = p0

! =========================
! EULER
! =========================

call cpu_time(start)

do i=0,N-1
    p_euler(i+1) = p_euler(i) + rb*(1.0d0 - p_euler(i))
end do

call cpu_time(finish)
time_euler = finish - start

! =========================
! TAYLOR ORDEN 2
! =========================

call cpu_time(start)

do i=0,N-1
    p_taylor(i+1) = p_taylor(i) + 0.001998d0*(1.0d0 - p_taylor(i))
end do

call cpu_time(finish)
time_taylor = finish - start

! =========================
! TRAPECIO
! =========================

call cpu_time(start)

do i=0,N-1
    p_trap(i+1) = (0.999d0*p_trap(i) + 0.002d0)/1.001d0
end do

call cpu_time(finish)
time_trap = finish - start

! =========================
! ERRORES ABSOLUTOS EN t=50
! =========================

error_euler = abs(p_exact(N) - p_euler(N))
error_taylor = abs(p_exact(N) - p_taylor(N))
error_trap = abs(p_exact(N) - p_trap(N))

! =========================
! GUARDAR ARCHIVOS .dat
! =========================

open(10, file="exacta.dat")
open(11, file="euler.dat")
open(12, file="taylor.dat")
open(13, file="trapecio.dat")

do i=0,N
    write(10,*) t(i), p_exact(i)
    write(11,*) t(i), p_euler(i)
    write(12,*) t(i), p_taylor(i)
    write(13,*) t(i), p_trap(i)
end do

close(10)
close(11)
close(12)
close(13)

! =========================
! RESULTADOS EN CONSOLA
! =========================

print *, "========== RESULTADOS =========="

print *
print *, "--- Euler ---"
print *, "Valor en t=50:", p_euler(N)
print *, "Error absoluto:", error_euler
print *, "Tiempo ejecucion:", time_euler, "segundos"

print *
print *, "--- Taylor Orden 2 ---"
print *, "Valor en t=50:", p_taylor(N)
print *, "Error absoluto:", error_taylor
print *, "Tiempo ejecucion:", time_taylor, "segundos"

print *
print *, "--- Trapecio ---"
print *, "Valor en t=50:", p_trap(N)
print *, "Error absoluto:", error_trap
print *, "Tiempo ejecucion:", time_trap, "segundos"

print *
print *, "Archivos generados: exacta.dat, euler.dat, taylor.dat, trapecio.dat"

end program metodos_modelo