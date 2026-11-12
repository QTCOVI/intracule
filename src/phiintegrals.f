C
C     Computes the integrals I(nx,ny) = integral_0^(2pi) sin^nx cos^ny'
C     for nx from 0 to NXMAX and ny from 0 to NYMAX
C
      SUBROUTINE phiintegrals (NXMAX, NYMAX, I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     Subrutina para calcular I(nx,ny) = integral_0^(2pi) sin^nx(phi) cos^ny(phi) dphi
C     con nx = 0..NXMAX, ny = 0..NYMAX
C
C     Argumentos:
C       NXMAX - (entrada) valor m�ximo de nx
C       NYMAX - (entrada) valor m�ximo de ny
C       I     - (salida)  arreglo bidimensional (0:NXMAX, 0:NYMAX) con los resultados
C
C     Nota: Se requiere que el llamante haya dimensionado I como (0:NXMAX, 0:NYMAX)
C
      DIMENSION I(0:NXMAX, 0:NYMAX)
      DIMENSION J(0:NXMAX, 0:NYMAX)
      DOUBLE PRECISION I, J, PI
C
C     Calcular pi
      PI = 4.0D0 * DATAN(1.0D0)
C
C     Inicializar matriz J (integral entre 0 y pi/2)
      DO NX = 0, NXMAX
         DO NY = 0, NYMAX
            J(NX, NY) = 0.0D0
         END DO
      END DO
C
C     Casos base para J(nx,ny) = int_0^(pi/2) sin^nx cos^ny
      J(0,0) = PI / 2.0D0
      IF (NXMAX .GE. 1) J(1,0) = 1.0D0
      IF (NYMAX .GE. 1) J(0,1) = 1.0D0
C
C     Rellenar usando recurrencia en nx (potencia del seno)
C     Se necesita J(nx-2,ny) para nx>=2
      DO NY = 0, NYMAX
         DO NX = 2, NXMAX
            J(NX, NY) = (DFLOAT(NX-1) / DFLOAT(NX+NY)) * J(NX-2, NY)
         END DO
      END DO
C
C     Rellenar usando recurrencia en ny (potencia del coseno)
C     Se necesita J(nx,ny-2) para ny>=2
      DO NX = 0, NXMAX
         DO NY = 2, NYMAX
            J(NX, NY) = (DFLOAT(NY-1) / DFLOAT(NX+NY)) * J(NX, NY-2)
         END DO
      END DO
C
C     Construir I(nx,ny) = integral_0^(2pi) sin^nx cos^ny
C     Solo es no nula si nx y ny son ambos pares
C     En ese caso: I = 4 * J(nx, ny)
      DO NX = 0, NXMAX
         DO NY = 0, NYMAX
            IF ((MOD(NX,2) .EQ. 0) .AND. (MOD(NY,2) .EQ. 0)) THEN
               I(NX, NY) = 4.0D0 * J(NX, NY)
            ELSE
               I(NX, NY) = 0.0D0
            ENDIF
         END DO
      END DO
C
      RETURN
      END
