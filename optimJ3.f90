program equation_solver
    implicit none
    
    real*8 :: A, B, C, f1, f2, f3, df1da, df1db, df1dc, df2da, df2db, df2dc, df3da, df3db, df3dc, delta
    integer :: max_iter, iter
    real*8, parameter :: epsilon = 1e-10
    
    ! Initial guesses for A, B, and C
    A = 2.660d0
    B = 2.2220d0
    C = 2.24d0
    
    ! Maximum number of iterations and initial iteration count
    max_iter = 100000
    iter = 0
    
    ! Iterate until convergence or maximum iterations reached
    do while (iter < max_iter)
        ! Evaluate the functions and their derivatives at the current values of A, B, and C
        f1 = 2.0d0*A + 5.0d0*B + 5.0d0*C - 2.0d0*sqrt(4.0d0*(B - C)**2 + (A - B)*(A - C)) - 0.697d0
        f2 = 5.0d0*A + 2.0d0*B + 5.0d0*C - 2.0d0*sqrt(4.0d0*(A - C)**2 - (A - B)*(B - C)) - 1.045d0
        f3 = 5.0d0*A + 5.0d0*B + 2.0d0*C - 2.0d0*sqrt(4.0d0*(A - B)**2 + (A - C)*(B - C)) - 1.094d0
        
        df1da = 2.0d0 - (2.0d0*(A - B)*(A - C)) / sqrt(4.0d0*(B - C)**2 + (A - B)*(A - C))
        df1db = 5.0d0 + (8.0d0*(B - C)) / sqrt(4.0d0*(B - C)**2 + (A - B)*(A - C))
        df1dc = 5.0d0 + (8.0d0*(B - C)) / sqrt(4.0d0*(B - C)**2 + (A - B)*(A - C))
        df2da = 5.0d0 + (8.0d0*(A - C)) / sqrt(4.0d0*(A - C)**2 - (A - B)*(B - C))
        df2db = 2.0d0 - (2.0d0*(A - C)*(B - C)) / sqrt(4.0d0*(A - C)**2 - (A - B)*(B - C))
        df2dc = 5.0d0 + (8.0d0*(A - C)) / sqrt(4.0d0*(A - C)**2 - (A - B)*(B - C))
        df3da = 5.0d0 + (8.0d0*(A - B)) / sqrt(4.0d0*(A - B)**2 + (A - C)*(B - C))
        df3db = 5.0d0 + (8.0d0*(A - B)) / sqrt(4.0d0*(A - B)**2 + (A - C)*(B - C))
        df3dc = 2.0d0 - (2.0d0*(A - C)*(B - C)) / sqrt(4.0d0*(A - B)**2 + (A - C)*(B - C))
        
        ! Calculate the determinant of the Jacobian matrix
        delta = df1da * (df2db * df3dc - df2dc * df3db) + df1db * (df2dc * df3da - df2da * df3dc)&
                + df1dc * (df2da * df3db - df2db * df3da)
        
        ! Update the values of A, B, and C using the Newton-Raphson method
        A = A - (f1 * (df2db * df3dc - df2dc * df3db) + f2 * (df1dc * df3db - df1db * df3dc)&
                + f3 * (df1db * df2dc - df1dc * df2db)) / delta
        B = B - (f1 * (df2dc * df3da - df2da * df3dc) + f2 * (df1da * df3dc - df1dc * df3da)&
                + f3 * (df1dc * df2da - df1da * df2dc)) / delta
        C = C - (f1 * (df2da * df3db - df2db * df3da) + f2 * (df1db * df3da - df1da * df3db)&
                + f3 * (df1da * df2db - df1db * df2da)) / delta
        open(10, file='J3AYDA', status='unknown')
        ! Print the values of the functions at each iteration
        write(*,*) "Iteration:", iter
        write(*,*) "f1 =", f1
        write(*,*) "f2 =", f2
        write(*,*) "f3 =", f3
        write(*,*) "------------------"
        
        ! Check for convergence
        if (abs(f1) < epsilon .and. abs(f2) < epsilon .and. abs(f3) < epsilon) then
            exit
        end if
        
        ! Increment the iteration count
        iter = iter + 1
    end do
    
    ! Print the results
    write(10,*) "Solution:"
    write(10,*) "A =", A,'cm-1'
    write(10,*) "B =", B,'cm-1'
    write(10,*) "C =", C,'cm-1'
    
end program equation_solver

