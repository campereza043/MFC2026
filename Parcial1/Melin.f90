PROGRAM MELIN
  IMPLICIT NONE

  ! ── Parámetros ────────────────────────
  INTEGER, PARAMETER :: N   = 10
  INTEGER, PARAMETER :: NPT = 600
  REAL(8), PARAMETER :: H    = 0.050D0
  REAL(8), PARAMETER :: ALFA = 0.2D0
  REAL(8), PARAMETER :: ANCHO = 1.0D0

  ! ── Variables ─────────────────────────
  REAL(8) :: Y(N), DYDT(N)
  REAL(8) :: T, DX
  REAL(8) :: t_inicio, t_fin, t_ejecucion
  INTEGER :: I

  DX = ANCHO / N
  T  = 0.0D0

  OPEN(UNIT=2, FILE='MELIN.dat', STATUS='REPLACE')

  ! ── Condición inicial ──────────────────
  DO I = 2, N
    Y(I) = 100.0D0
  END DO

  ! ── Condición de frontera en X=0 ───────
  Y(1) = 0.0D0

  ! ── Cronómetro inicia ──────────────────
  CALL CPU_TIME(t_inicio)

  CALL SOL(T, Y, DYDT, H, N, NPT, ALFA, DX)

  ! ── Cronómetro termina ─────────────────
  CALL CPU_TIME(t_fin)
  t_ejecucion = t_fin - t_inicio

  CLOSE(2)

  WRITE(*,'(A)')         '================================='
  WRITE(*,'(A,F10.6,A)') '  Tiempo de ejecucion: ', t_ejecucion, ' segundos'
  WRITE(*,'(A)')         '================================='

  STOP

CONTAINS

  ! =============================================
  SUBROUTINE DERIVS(X, Y, DYDT, ALFA, DX, N)
    IMPLICIT NONE
    INTEGER,  INTENT(IN)  :: N
    REAL(8),  INTENT(IN)  :: X, ALFA, DX
    REAL(8),  INTENT(IN)  :: Y(N)
    REAL(8),  INTENT(OUT) :: DYDT(N)

    REAL(8) :: ALFA2, DX2, FR_DER
    INTEGER :: I

    ALFA2   = ALFA * ALFA
    DX2     = DX * DX
    DYDT(1) = 0.0D0

    ! Puntos interiores
    DO I = 2, N-1
      DYDT(I) = ALFA2 * (Y(I+1) - 2.0D0*Y(I) + Y(I-1)) / DX2
    END DO

    ! Frontera aislada x=l  →  Neumann dT/dx=0
    FR_DER  = (4.0D0*Y(N) - Y(N-1)) / 3.0D0
    DYDT(N) = ALFA2 * (FR_DER - 2.0D0*Y(N) + Y(N-1)) / DX2

  END SUBROUTINE DERIVS

  ! =============================================
  SUBROUTINE RK4(Y, DYDT, N, X, H, YOUT, ALFA, DX)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: N
    REAL(8), INTENT(IN)  :: Y(N), DYDT(N)
    REAL(8), INTENT(IN)  :: X, H, ALFA, DX
    REAL(8), INTENT(OUT) :: YOUT(N)

    REAL(8) :: YT(N), DYT(N), DYM(N)
    REAL(8) :: HH, H6, XH
    INTEGER :: I

    HH = H * 0.5D0
    H6 = H / 6.0D0
    XH = X + HH

    ! ── Paso 1 ────────────────────────────
    DO I = 1, N
      YT(I) = Y(I) + HH * DYDT(I)
    END DO
    CALL DERIVS(XH, YT, DYT, ALFA, DX, N)

    ! ── Paso 2 ────────────────────────────
    DO I = 1, N
      YT(I) = Y(I) + HH * DYT(I)
    END DO
    CALL DERIVS(XH, YT, DYM, ALFA, DX, N)

    ! ── Paso 3 ────────────────────────────
    DO I = 1, N
      YT(I)  = Y(I) + H * DYM(I)
      DYM(I) = DYT(I) + DYM(I)
    END DO
    CALL DERIVS(X+H, YT, DYT, ALFA, DX, N)

    ! ── Combinación final RK4 ─────────────
    DO I = 1, N
      YOUT(I) = Y(I) + H6 * (DYDT(I) + DYT(I) + 2.0D0*DYM(I))
    END DO

  END SUBROUTINE RK4

  ! =============================================
  SUBROUTINE SOL(T, Y, DYDT, H, N, NPT, ALFA, DX)
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: N, NPT
    REAL(8), INTENT(INOUT) :: T, Y(N), DYDT(N)
    REAL(8), INTENT(IN)    :: H, ALFA, DX

    REAL(8) :: YOUT(N), FR_DER
    INTEGER :: I, J

    ! ── Primer paso ───────────────────────
    FR_DER = (4.0D0*Y(N) - Y(N-1)) / 3.0D0
    WRITE(2,'(15F10.4)') T, (Y(I), I=1,N), FR_DER

    CALL DERIVS(T, Y, DYDT, ALFA, DX, N)
    CALL RK4(Y, DYDT, N, T, H, YOUT, ALFA, DX)

    FR_DER = (4.0D0*YOUT(N) - YOUT(N-1)) / 3.0D0
    WRITE(2,'(15F10.4)') T+H, (YOUT(I), I=1,N), FR_DER
    T = T + H

    ! ── Ciclo principal ───────────────────
    DO J = 1, NPT
      T = T + H
      DO I = 1, N
        Y(I) = YOUT(I)
      END DO
      CALL DERIVS(T, Y, DYDT, ALFA, DX, N)
      CALL RK4(Y, DYDT, N, T, H, YOUT, ALFA, DX)
      FR_DER = (4.0D0*YOUT(N) - YOUT(N-1)) / 3.0D0
      WRITE(2,'(15F10.4)') T, (YOUT(I), I=1,N), FR_DER
    END DO

  END SUBROUTINE SOL

END PROGRAM MELIN