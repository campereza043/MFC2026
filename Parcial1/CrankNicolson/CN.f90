PROGRAM CRANK_NICOLSON
  IMPLICIT NONE

  ! ── Parámetros ────────────────────────────────────────────
  INTEGER, PARAMETER :: N    = 10
  INTEGER, PARAMETER :: NPT  = 600
  REAL(8), PARAMETER :: H    = 0.050D0
  REAL(8), PARAMETER :: ALFA = 0.2D0
  REAL(8), PARAMETER :: ANCHO = 1.0D0

  ! ── Variables ─────────────────────────────────────────────
  REAL(8) :: Y(N), Y_NEW(N)
  REAL(8) :: A(N), B(N), C(N), D(N)   ! sistema tridiagonal
  REAL(8) :: T, DX, R
  REAL(8) :: t_inicio, t_fin, t_ejecucion
  REAL(8) :: FR_DER
  INTEGER :: I, J

  DX = ANCHO / N
  T  = 0.0D0
  R  = ALFA*ALFA * H / (2.0D0 * DX*DX)   ! número de Fourier/2

  OPEN(UNIT=3, FILE='CNF.dat', STATUS='REPLACE')

  ! ── Condición inicial ──────────────────────────────────────
  DO I = 2, N
    Y(I) = 100.0D0
  END DO
  Y(1) = 0.0D0   ! frontera fría fija

  ! ── Cronómetro ────────────────────────────────────────────
  CALL CPU_TIME(t_inicio)

  ! ── Escribir t=0 ──────────────────────────────────────────
  FR_DER = (4.0D0*Y(N) - Y(N-1)) / 3.0D0
  WRITE(3,'(15F10.4)') T, (Y(I), I=1,N), FR_DER

  ! ── Ciclo temporal ────────────────────────────────────────
  DO J = 1, NPT

    ! Construir sistema tridiagonal A·Y_NEW = D
    ! Frontera x=0: Dirichlet T1=0 siempre
    A(1) = 0.0D0
    B(1) = 1.0D0
    C(1) = 0.0D0
    D(1) = 0.0D0

    ! Puntos interiores i=2..N-1
    DO I = 2, N-1
      A(I) = -R
      B(I) =  1.0D0 + 2.0D0*R
      C(I) = -R
      D(I) =  R*Y(I-1) + (1.0D0 - 2.0D0*R)*Y(I) + R*Y(I+1)
    END DO

    ! Frontera x=l: Neumann dT/dx=0
    ! T_{N+1} = (4*T_N - T_{N-1})/3  →  se sustituye en la ecuación
    ! La fila N queda:
    !   -(R + 4R/3)*Y_NEW(N-1) + ... no, se usa la aproximación directa:
    ! Con T_{N+1} = (4T_N - T_{N-1})/3 la ecuación para el nodo N es:
    ! dT_N/dt = alfa²*(T_{N+1} - 2T_N + T_{N-1})/dx²
    !         = alfa²*((4T_N-T_{N-1})/3 - 2T_N + T_{N-1})/dx²
    !         = alfa²*(T_{N-1}/3 - 2T_N/3)/dx²
    ! Aplicando Crank-Nicolson:
    ! (Y_NEW(N)-Y(N))/H = (alfa²/dx²)*0.5*(
    !   (Y(N-1)/3 - 2Y(N)/3) + (Y_NEW(N-1)/3 - 2Y_NEW(N)/3) )
    ! Reordenando:
    !   -R/3 * Y_NEW(N-1)  +  (1 + R*2/3) * Y_NEW(N)  =
    !    R/3 * Y(N-1)      +  (1 - R*2/3) * Y(N)
    A(N) = -R / 3.0D0
    B(N) =  1.0D0 + (2.0D0*R) / 3.0D0
    C(N) =  0.0D0
    D(N) =  (R/3.0D0)*Y(N-1) + (1.0D0 - (2.0D0*R)/3.0D0)*Y(N)

    ! ── Eliminación de Thomas (TDMA) ──────────────────────
    DO I = 2, N
      B(I) = B(I) - A(I) * C(I-1) / B(I-1)
      D(I) = D(I) - A(I) * D(I-1) / B(I-1)
    END DO

    ! Sustitución regresiva
    Y_NEW(N) = D(N) / B(N)
    DO I = N-1, 1, -1
      Y_NEW(I) = (D(I) - C(I)*Y_NEW(I+1)) / B(I)
    END DO

    T = T + H
    DO I = 1, N
      Y(I) = Y_NEW(I)
    END DO

    FR_DER = (4.0D0*Y(N) - Y(N-1)) / 3.0D0
    WRITE(3,'(15F10.4)') T, (Y(I), I=1,N), FR_DER

  END DO

  CALL CPU_TIME(t_fin)
  t_ejecucion = t_fin - t_inicio

  CLOSE(3)

  WRITE(*,'(A)')         '================================='
  WRITE(*,'(A)')         '  MELIN — Crank-Nicolson Fortran'
  WRITE(*,'(A)')         '---------------------------------'
  WRITE(*,'(A,F8.4)')    '  R (Fourier/2) = ', R
  WRITE(*,'(A)')         '  Estabilidad   : INCONDICIONAL'
  WRITE(*,'(A,F10.6,A)') '  Tiempo        : ', t_ejecucion, ' s'
  WRITE(*,'(A)')         '  Archivo       : MELINCN.dat'
  WRITE(*,'(A)')         '================================='

  STOP
END PROGRAM CRANK_NICOLSON