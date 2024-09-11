module markovSparse
    use fsparse
    use denseMatrix
    use iso_fortran_env, only: dp=>real64
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
            real(dp), dimension(SIZE(probabilityCoeffs)) :: prevProbCoeffs, timeStepUntilConvergenceSparse
            real(dp) :: tol
            logical converged
            integer :: timeSteps, stepsTaken
            
            stepsTaken = 0
            tol = 1.0/(SIZE(probabilityCoeffs) * 1e6) !Tolerance is one part in one million.
            timeSteps = 20
            tol = tol * timeSteps !Dynamically change tolerance based on number of steps taken
            
            converged = .FALSE.
            prevProbCoeffs = probabilityCoeffs
            do while (.not. converged)
                timeStepUntilConvergenceSparse = timeStepSparse(matrix, prevProbCoeffs, timeSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
                converged = convergence(timeStepUntilConvergenceSparse, prevProbCoeffs, tol)
                prevProbCoeffs = timeStepUntilConvergenceSparse
                stepsTaken = stepsTaken + timeSteps
                if(mod(stepsTaken, 500) == 0) then
                    print *, 'Steps taken: ', stepsTaken
                end if
            end do
            print*, 'Steps taken: ', stepsTaken
        end function timeStepUntilConvergenceSparse


        function timeStep(matrix, probabilityCoeffs, steps) !This computes repeated matrix multiplication (A*(A*...*(A*V)
            integer, intent(in) :: steps                     !In order to only have to store one matrix
            double precision, dimension(:), intent(in) :: probabilityCoeffs
            double precision, dimension(SIZE(probabilityCoeffs)) :: timeStep
            double precision, dimension(SIZE(probabilityCoeffs), SIZE(probabilityCoeffs)), intent(in):: matrix 
            integer :: i
            timeStep = probabilityCoeffs
            do i = 1, steps
                timeStep = MATMUL(matrix, timeStep)
            end do
            
        end function timeStep

        function timeStepUntilConvergence(matrix, probabilityCoeffs)
            real(dp), dimension(:), intent(in) :: probabilityCoeffs
            real(dp), dimension(SIZE(probabilityCoeffs),SIZE(probabilityCoeffs)), intent(in):: matrix
            real(dp), dimension(SIZE(probabilityCoeffs)) :: prevProbCoeffs, timeStepUntilConvergence
            real(dp) :: tol
            logical converged
            integer :: timeSteps, stepsTaken
            
            stepsTaken = 0
            tol = 1.0/(SIZE(probabilityCoeffs) * 1e6) !Tolerance is one part in one million.
            timeSteps = 20
            tol = tol * timeSteps !Dynamically change tolerance based on number of steps taken
            
            converged = .FALSE.
            prevProbCoeffs = probabilityCoeffs
            do while (.not. converged)
                timeStepUntilConvergence = timeStep(matrix, prevProbCoeffs, timeSteps) !Dont check for convergence after every step, rather take a few time steps at a time.
                converged = convergence(timeStepUntilConvergence, prevProbCoeffs, tol)
                prevProbCoeffs = timeStepUntilConvergence
                stepsTaken = stepsTaken + timeSteps
                if(mod(stepsTaken, 500) == 0) then
                    print *, 'Steps taken: ', stepsTaken
                end if
            end do
            print*, 'Steps taken: ', stepsTaken
            
        end function timeStepUntilConvergence
        
        !TODO
        
        function convergence(newCoeff, oldCoeff, tol)
            !Checks if all elements in two same sized vectors are close neough given a tolerance
            real(dp), dimension(:), intent(in) :: newCoeff, oldCoeff
            real(dp), dimension(SIZE(newCoeff)) :: difference
            real(dp), intent(in) :: tol
            LOGICAL :: convergence
            integer i
            convergence = .TRUE.
            difference = abs(newCoeff - oldCoeff)
            do i = 1, SIZE(difference)
                if(difference(i) > tol) then
                    convergence = .FALSE.
                    return
                end if
            end do
        end function convergence
end module markovSparse