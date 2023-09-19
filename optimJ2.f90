program equation_solver
    implicit none
    
    real(8) :: A, B, C, f1, f2, f3, df1da, df1db, df1dc, df2da, df2db, df2dc, df3da, df3db, df3dc, delta
    integer :: max_iter, iter
    real(8), parameter :: epsilon = 1e-9
    
    ! Initial guesses for A, B, and C (updated)
    A = 3.3345d0
    B = 2.2229d0
    C = 2.2226d0
    
    ! Maximum number of iterations and initial iteration count
    max_iter = 1000
    iter = 0
    
    ! Iterate until convergence or maximum iterations reached
    do while (iter < max_iter)
        ! Evaluate the functions and their derivatives at the current values of A, B, and C
        f1 = A + B + 4.0d0*C - 0.708d0
        f2 = A + 4.0*B + C - 0.733d0
        f3 = 2.0d0*A + 2.0d0*B + 2.0d0*C - 2.0d0*dsqrt((B - C)**2 + (A - C)*(A - B)) - 0.349d0
        df1da = 1.0d0
        df1db = 1.0d0
        df1dc = 4.0d0
        df2da = 1.0d0
        df2db = 4.0d0
        df2dc = 1.0d0
        df3da = 2.0d0 - (2.0d0*(A - C)) / sqrt((B - C)**2 + (A - C)*(A - B))
        df3db = 2.0d0 + (B - C) / sqrt((B - C)**2 + (A - C)*(A - B))
        df3dc = 2.0d0 - (2.0d0*(C - A)) / sqrt((B - C)**2 + (A - C)*(A - B))
        
        ! Calculate the determinant of the Jacobian matrix
        delta = df1da * (df2db * df3dc - df2dc * df3db) + df1db &
        * (df2dc * df3da - df2da * df3dc) + df1dc * (df2da * df3db - df2db * df3da)
        
        ! Update the values of A, B, and C using the Newton-Raphson method
        A = A - (f1 * (df2db * df3dc - df2dc * df3db) &
        + f2 * (df1dc * df3db - df1db * df3dc) + f3 * (df1db * df2dc - df1dc * df2db)) / delta
        B = B - (f1 * (df2dc * df3da - df2da * df3dc) &
        + f2 * (df1da * df3dc - df1dc * df3da) + f3 * (df1dc * df2da - df1da * df2dc)) / delta
        C = C - (f1 * (df2da * df3db - df2db * df3da)&
         + f2 * (df1db * df3da - df1da * df3db) + f3 * (df1da * df2db - df1db * df2da)) / delta
       open(10, file='J2Ayda', status='unknown')
        ! Output current iteration information
        write(*, *) "Iteration: ", iter
        write(*, *) "f1 = ", f1
        write(*, *) "f2 = ", f2
        write(*, *) "f3 = ", f3
        write(*, *) "------------------"
        
        ! Check for convergence
        if (abs(f1) < epsilon .and. abs(f2) < epsilon .and. abs(f3) < epsilon) then
            exit
        end if
        
        ! Increment the iteration count
        iter = iter + 1
    end do
    
    ! Print the results
    write(10, *) "Solution:"
    write(10, *) "A = ", A, " cm-1"
    write(10, *) "B = ", B, " cm-1"
    write(10, *) "C = ", C, " cm-1"
    
end program equation_solver

