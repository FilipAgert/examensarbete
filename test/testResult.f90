program testResult
    use iso_fortran_env, only: dp=>real64
    use result_module
    implicit none

    real(dp),dimension (5,2) :: probabilityDensity
    real, dimension(2) :: solveTime
    integer, dimension(2) :: matrixMultiplications
    real, dimension(2) :: fusionFraction
    real, dimension(2) :: fissionFraction
    real, dimension(3,2) :: startCoordinate
    real, dimension(2) :: energies
    type(convergedResult) :: res
    ! Initialize arrays
    probabilityDensity = reshape( (/ 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
                                     0.1_dp, 0.5_dp, 0.2_dp, 0.4_dp, 0.5_dp /), &
                                  shape(probabilityDensity) )

    solveTime = (/ 10.0, 12.0 /)
    matrixMultiplications = (/ 3, 5 /)
    fusionFraction = (/ 0.6, 0.7 /)
    fissionFraction = (/ 0.4, 0.3 /)
    startCoordinate = reshape( (/ 1.0, 2.0, 3.0, 4.0, 5.0, &
                                 5.0, 4.0, 3.0, 2.0, 1.0 /), &
                               shape(startCoordinate) )
    energies = (/ 1.0, 2.0/)

    call res%addResult(probabilityDensity(:,1), startCoordinate(:,1), energies(1), solveTime(1),&
                matrixMultiplications(1), fusionFraction(1), fissionFraction(1))
    call res%addResult(probabilityDensity(:,2), startCoordinate(:,2), energies (2), solveTime(2),& 
                matrixMultiplications(2), fusionFraction(2), fissionFraction(2))

    call res%printResult()
    

end program testResult