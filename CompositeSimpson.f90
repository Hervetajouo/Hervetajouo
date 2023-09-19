REAL(8) FUNCTION f(x) RESULT(Func)
      REAL(8), INTENT(IN) :: x 
      Func =x*LOG(x)
      END FUNCTION f 

      PROGRAM CompSimpson
              IMPLICIT NONE 
              REAL(8) ::x, a, b, h, S0,  S1, S2, Integ 
              REAL(8), EXTERNAL :: f 
              INTEGER :: n, i 
               a = 1.0 ; b = 2.0 ; n = 4
               h = (b - a )/n
               S0 = f(a) + f(b) 
               S1 = 0 ; S2= 0 
               DO i = 1,( n-1) 
               x = a + i*h  
               IF ( MOD(i, 2)==0 ) THEN 
                       S1 = S1 + f(x)
               ELSE
                       S2 = S2 + f(x)
               END IF 
               END DO 
               Integ = h*( S0 + 2*S1 + 4*S2 )/3.0
               PRINT*, 'The value of the integral is', Integ 
               END PROGRAM CompSimpson

