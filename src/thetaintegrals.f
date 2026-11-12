C     PROGRAM PRUEBA_INTEGRALES
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     Programa de prueba para la subrutina CALCULA_INTEGRALES
C
C     PARAMETER (maxnxy = 20, maxnz = 20)
C     DIMENSION XI(0:maxnxy, 0:maxnz)
C     DOUBLE PRECISION XI
C     INTEGER P,Q
C
C     Llamar a la subrutina
C     CALL thetaintegrals (maxnxy, maxnz, XI)
C
C     Escribir resultados
C     WRITE(6,*) 'Resultados I(nxy,nz) = int0^pi sin^(nxy+1) cos^nz'
C     WRITE(6,*) 'nxy filas, nz columnas'
C     WRITE(6,*) '----------------------------------------'
C     WRITE(6,100) (NZ, NZ=0,maxnz)
C 100 FORMAT(11X, 11(I3,7X))
C
C     DO 10 NXY = 0, maxnxy
C        WRITE(6,110) NXY, (XI(NXY, NZ), NZ=0,maxnz)
C 110    FORMAT('nxy=',I2, ' | ', 29(F8.5,1X))
C  10 CONTINUE
C
C     STOP
C     END

C     Computes the integrals I(nxy,nz) = int_0^pi sin^(nxy+1) cos^nz
C     for nxy=0,...,maxnxy and nz=0,...maxnz

      SUBROUTINE thetaintegrals (maxnxy, maxnz, XI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     Subrutina para calcular 
C     I(nxy,nz) = integral_0^pi sin^(nxy+1)(theta) cos^nz(theta) dtheta
C     con nxy = 0..maxnxy, nz = 0..maxnz
C
C     Argumentos:
C       maxnxy  - (entrada) valor m�ximo de nxy
C       maxnz  - (entrada) valor m�ximo de nz
C       I  - (salida)  arreglo bidimensional (0:maxnxy, 0:maxnz) con los resultados
C
C     Nota: Se requiere que el llamante haya dimensionado I como
C     (0:maxnxy, 0:maxnz)
C
      DIMENSION XI(0:maxnxy, 0:maxnz)
      DIMENSION XJ(0:maxnxy+1, 0:maxnz)
      DOUBLE PRECISION XI, XJ, PI
      INTEGER P,Q
C
C     Calcular pi
      PI = 4.0D0 * DATAN(1.0D0)
C
C     Inicializar matriz J (integral entre 0 y pi/2)
      DO P = 0, maxnxy+1
         DO Q = 0, maxnz
            XJ(P, Q) = 0.0D0
         END DO
      END DO
C
C     Casos base para J(p,q) = int_0^(pi/2) sin^p cos^q
      XJ(0,0) = PI / 2.0D0
      IF (maxnxy+1 .GE. 1) XJ(1,0) = 1.0D0
      IF (maxnz .GE. 1) XJ(0,1) = 1.0D0
C
C     Rellenar usando recurrencia en p (potencia del seno)
C     Se necesita J(p-2,q) para p>=2
      DO Q = 0, maxnz
         DO P = 2, maxnxy+1
            XJ(P, Q) = (DFLOAT(P-1) / DFLOAT(P+Q)) * XJ(P-2, Q)
         END DO
      END DO
C
C     Rellenar usando recurrencia en q (potencia del coseno)
C     Se necesita J(p,q-2) para q>=2
      DO P = 0, maxnxy+1
         DO Q = 2, maxnz
            XJ(P, Q) = (DFLOAT(Q-1) / DFLOAT(P+Q)) * XJ(P, Q-2)
         END DO
      END DO
C
C     Construir I(nxy,nz) = integral_0^pi sin^(nxy+1) cos^nz
C     Si nz es par -> I = 2*J(nxy+1, nz); si nz impar -> I = 0
      DO NXY = 0, maxnxy
         DO NZ = 0, maxnz
            P = NXY + 1
            IF (MOD(NZ,2) .EQ. 0) THEN
               XI(NXY, NZ) = 2.0D0 * XJ(P, NZ)
            ELSE
               XI(NXY, NZ) = 0.0D0
            ENDIF
*           write (66,*) nxy,nz,xi(nxy,nz)
         END DO
      END DO
C
      RETURN
      END
