REAL(8) FUNCTION f(x) RESULT(Func)
      REAL(8), INTENT(IN) :: x
      Func =x*LOG(x)
      END FUNCTION f

      PROGRAM CompTrapez
              IMPLICIT NONE
              REAL(8) ::x, a, b, h, S0,  S1, Integ
              REAL(8), EXTERNAL :: f
              INTEGER :: n, i
               a = 1.0 ; b = 2.0 ; n = 4
               h = (b - a )/n
               S0 = f(a) + f(b)
               DO i = 1, (n-1)
               x = a + i*h
               S1 = S1 + f(x)
               END DO
               Integ = h*(S0 + 2*S1)/2.0
               PRINT*, ' The value of the Integral is,' , Integ 
               END PROGRAM CompTrapez

