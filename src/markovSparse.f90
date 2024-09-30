module markovSparse
    use fsparse
    use denseMatrix
    use iso_fortran_env, only: dp=>real64, dp=>real32
    implicit none

    contains
        

        function timeStepSparse(matrixSparse, probabilityCoeffs, steps)
            integer, intent(in) :: steps
            real(dp), dimension(:), intent(in) :: probabilityCoeffs
            real(dp), dimension(SIZE(probabilityCoeffs)) :: timeStepSparse, temp
            type(COO_dp) :: matrixSparse
            integer :: i
            timeStepSparse = probabilityCoeffs
            do i = 1, steps
                temp = 0
                call matvec(matrixSparse, timeStepSparse, temp)
                timeStepSparse = temp ! new timeStepSparse = matrixSparse * timeStepSparse
            end do
        end function timeStepSparse
        
        
        function timeStepUntilConvergenceSparse(matrix, probabilityCoeffs)
            real(dp), dimension(:), intent(in) :: probabilityCoeffs
            type(COO_dp) :: matrix
            real(dp), dimension(SIZE(probabilityCoeffs)) :: prevProbCoeffs, timeStepUntilConvergenceSparse, diff
            real(dp) :: tol, normVal
            logical converged
            integer :: timeSteps, matrixMultiplications
            character (len = 2):: normType = "L2"
            
            matrixMultiplications = 0
            tol = 1.0/(SIZE(probabilityCoeffs) * 1e6) !Tolerance is one part in one million.
            timeSteps = 20
            tol = tol * timeSteps !Dynamically change tolerance based on number of matrix multiplications in a row
            
            converged = .FALSE.
            prevProbCoeffs = probabilityCoeffs
            do while (.not. converged)
                timeStepUntilConvergenceSparse = timeStepSparse(matrix, prevProbCoeffs, timeSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
                converged = convergence(timeStepUntilConvergenceSparse, prevProbCoeffs, tol, normType)
               
                if(mod(matrixMultiplications, 500) == 0) then

                    diff = timeStepUntilConvergenceSparse - prevProbCoeffs
                    normVal = norm(diff, normType)
                    print *, 'Matrix multiplications: ', matrixMultiplications, " ", normType, " norm: ", normVal
                end if
                prevProbCoeffs = timeStepUntilConvergenceSparse
                matrixMultiplications = matrixMultiplications + timeSteps
            end do
            print*, 'Matrix multiplications: ', matrixMultiplications
        end function timeStepUntilConvergenceSparse


        function timeStep(matrix, probabilityCoeffs, steps) !This computes repeated matrix multiplication (A*(A*...*(A*V)
            integer, intent(in) :: steps                     !In order to only have to store one matrix
            real(dp), dimension(:), intent(in) :: probabilityCoeffs
            real(dp), dimension(SIZE(probabilityCoeffs)) :: timeStep
            real(dp), dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in):: matrix 
            integer :: i
            timeStep = probabilityCoeffs
            do i = 1, steps
                timeStep = MATMUL(matrix, timeStep)
            end do
            
        end function timeStep

        function timeStepUntilConvergence(matrix, probabilityCoeffs)
            real(dp), dimension(:), intent(in) :: probabilityCoeffs
            real(dp), dimension(SIZE(probabilityCoeffs),SIZE(probabilityCoeffs)), intent(in):: matrix
            real(dp), dimension(SIZE(probabilityCoeffs)) :: prevProbCoeffs, timeStepUntilConvergence, diff
            real(dp) :: tol, normVal
            logical converged
            integer :: timeSteps, matrixMultiplications
            character (len = 2):: normType = "L2"

            matrixMultiplications = 0
            tol = 1.0/(SIZE(probabilityCoeffs) * 1e6) !Tolerance is one part in one million.
            timeSteps = 20
            tol = tol * timeSteps !Dynamically change tolerance based on number of steps taken
            
            converged = .FALSE.
            prevProbCoeffs = probabilityCoeffs
            do while (.not. converged)
                timeStepUntilConvergence = timeStep(matrix, prevProbCoeffs, timeSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
                converged = convergence(timeStepUntilConvergence, prevProbCoeffs, tol, normType)
                prevProbCoeffs = timeStepUntilConvergence
                matrixMultiplications = matrixMultiplications + timeSteps
                if(mod(matrixMultiplications, 500) == 0) then
                    diff = timeStepUntilConvergence - prevProbCoeffs
                    normVal = norm(diff, normType)
                    print *, 'Matrix multiplications: ', matrixMultiplications, " ", normType, " norm: ", normVal
                    
                end if
            end do
            print*, 'Matrix multiplications: ', matrixMultiplications
            
        end function timeStepUntilConvergence
        
        !TODO
        
        function convergence(newCoeff, oldCoeff, tol, normType)
            !Checks if all elements in two same sized vectors are close neough given a tolerance
            real(dp), dimension(:), intent(in) :: newCoeff, oldCoeff
            real(dp), dimension(SIZE(newCoeff)) :: difference
            real(dp), intent(in) :: tol
            real(dp) :: normVal
            LOGICAL :: convergence
            character (len = *), intent(in), optional :: normType
            convergence = .TRUE.
            difference = (newCoeff - oldCoeff)
            normVal = norm(difference, normType)
            convergence = (normVal < tol)
        end function convergence

        function norm(vec, type)
            real(dp), dimension(:), intent(in) :: vec
            character (len = *), optional :: type !Set Linf if you want L_inf norm. Default is L2 norm
            real(dp) :: norm
            select case (type)
                case ('Linf')
                    norm = MAXVAL(abs(vec)) !L_infinity norm. Returns max value of vector
                case default !L_2 norm
                    norm = NORM2(vec)
            end select
        end function norm
end module markovSparse