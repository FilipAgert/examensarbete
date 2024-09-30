module linearInterpolation
    ! use stdlib_linalg, only: Operator(.det.), Operator(.inv.), eig
    use iso_fortran_env, only: dp=>real64, sp=>real32
    !use stdlib_kinds, only: sp
    use denseMatrix, only: printMatrixS, printMatrixToFile
    implicit none
    public guessPdfFast
contains
    ! function pdfInterpolation(coord, E, E0, E1, pdf0, pdf1) result (interpolation)
    !     real, intent(in) :: coord(:,:), E, E0, E1, pdf0(:), pdf1(:)
    !     real :: interpolation(size(pdf0))
    !     real :: alpha, maxDeltaMu, kDet
    !     real, dimension(size(coord,1)) :: mu0, mu1
    !     real, dimension(size(coord,1), size(coord,1)) :: K, K0, K1, muMult
    !     real, dimension(size(coord,1),1) :: deltaMu
    !     alpha = (E-E0)/(E1-E0)  !Some measure of how close the guess is to previous guesses.

    !     mu0 = getMean(coord, pdf0)
    !     mu1 = getMean(coord, pdf1)
    !     deltaMu(:,1) = abs(mu0-mu1)
    !     muMult = MATMUL(deltaMu, transpose(deltaMu))
    !     maxDeltaMu = maxval(muMult)

    !     K0 = getCovariance(coord, pdf0)
    !     K1 = getCovariance(coord, pdf1)
    !     K = (1-alpha)*K0 + alpha*K1
    !     kDet = .det.K
        

    !     interpolation = 0
    !     if(maxDeltaMu * 100 < kDet) then !if maxDeltaMu is much smaller than determinant then fast method is viable.
    !         interpolation = guessPdfFast(coord, E, E0, E1, pdf0, pdf1)
    !     else !else do slow method.
    !         interpolation = guessPdfSlow(coord, E, E0, E1, pdf0, pdf1)
    !     endif
        
    ! end function pdfInterpolation

    ! function sqrtm(A) result(sqrtA)
    !     real, intent(inout) :: A(:,:)
    !     real, allocatable :: sqrtA(:,:), D(:,:), lambda(:)
    !     real, allocatable :: S(:,:), Sinv(:,:)
    !     complex, allocatable :: Sc(:,:),lambdac(:)
    !     integer :: n, i
    !     n = size(A,1)
        
    !     allocate(lambdac(n), lambda(n), Sc(n,n), S(n,n), D(n,n), Sinv(n,n), sqrtA(n,n))
    !     D = 0
    !     call eig(A,lambdac, right=Sc)
    !     lambda = real(lambdac)
    !     S = real(Sc)
    !     do i = 1,SIZE(lambda)
    !         D(i,i) = sqrt(lambda(i))
    !     end do
    !     Sinv = .inv.S
    !     sqrtA = MATMUL(S,MATMUL(D,Sinv))
    ! end function sqrtm

    ! function guessPdfSlow(coord, E, E0, E1, pdf0, pdf1) result(interpolation)
    !     real, intent(in) :: coord(:,:), E, E0, E1, pdf0(:), pdf1(:)
    !     real(dp) :: interpolation(size(pdf0))
    !     real :: alpha, sumInt
    !     real, dimension(size(coord,1)) :: mu0, mu1, mu
    !     real, dimension(size(coord,1), size(coord,1)) :: K, K0, K1, muMult, K0Kinv, K1Kinv, Kinv
    !     integer, dimension(size(coord,2), size(coord,2)) :: gamma0, gamma1
    !     integer, dimension(size(coord,2)) :: s0, s1
    !     integer :: i
    !     real, dimension(size(coord,1), size(coord,2)) :: x0,x1
    !     alpha = (E-E0)/(E1-E0)
    !     mu0 = getMean(coord, pdf0)

    !     mu1 = getMean(coord, pdf1)
    !     mu = (1-alpha)*mu0 + alpha*mu1
    !     K0 = getCovariance(coord, pdf0)
    !     K1 = getCovariance(coord, pdf1)
    !     K = (1-alpha)*K0 + alpha*K1

    !     Kinv = .inv.K
    !     K0Kinv = MATMUL(K0, Kinv)
    !     K1Kinv = MATMUL(K1, Kinv)


    !     x0 = MATMUL(sqrtm(K0Kinv),(coord-spread(mu,2,SIZE(coord,2)))) + spread(mu0,2,SIZE(coord,2))
    !     x1 = MATMUL(sqrtm(K1Kinv),(coord-spread(mu,2,SIZE(coord,2)))) + spread(mu1,2,SIZE(coord,2))

    !     call generateGammaMatrix(gamma0, s0, coord, x0, 10)
    !     call generateGammaMatrix(gamma1, s1, coord, x1, 10)
    !     call printMatrixToFile("gamma", real(gamma0,8))
    !     interpolation = 0
    !     do i = 1,size(coord,2)
            
    !         interpolation = interpolation + (1-alpha)*real(gamma0(i,:))*pdf0(i)/real(s0(i))
    !         interpolation = interpolation + alpha*real(gamma1(i,:))*pdf1(i)/real(s1(i))
    !     end do
    !     sumInt = sum(interpolation)
    !     interpolation = interpolation/sumInt

    ! end function guessPdfSlow

    function guessPdfFast(coord, E,E0,E1, pdf0, pdf1) result(interpolation)
        real(sp), intent(in) :: coord(:,:), E, E0, E1
        real(dp), intent(in) :: pdf0(:), pdf1(:)
        real(dp) :: interpolation(size(pdf0))
        real(dp) :: alpha, sumInt
        alpha = (E-E0)/(E1-E0)
        interpolation = (1-alpha)*pdf0 + alpha*pdf1
        sumInt= sum(interpolation)
        interpolation = interpolation/sumInt
    end function guessPdfFast

!     function getMean(x,val) result(mean)
!         real, intent(in) :: x(:,:), val(:) 
!         !x is a set of coordinates of size DxN, where D is number of dimensions, and N is number of points
!         !Val is the associated  probability distribution of size 1xN
!         real :: mean(size(x,1)), coord(size(x,1)) !Return value is of size Dx1. Returns mean for each dimension of input.
!         integer :: i
!         mean = 0
!         do i = 1, size(x,2)
!             coord = x(:,i)
!             mean = mean + coord*val(i)
!         end do
!     end function getMean

!     function getCovariance(x,val) result(covariance)
!         real, intent(in) :: x(:,:), val(:) 
!         !Let val be 1-dimensional of length M
!         !Let x be DxM dimensional
!         !Where D is number of dimensions.
!         real :: covariance(size(x,1),size(x,1)), coord(size(x,1)), mean(size(x,1)), u(size(x,1),1)
!         integer :: i
!         mean = getMean(x, val)
!         covariance = 0
!         do i = 1,size(x,2)
!             coord = x(:,i)
!             u(:,1) = coord-mean
!             covariance = covariance + MATMUL(u, TRANSPOSE(u))*val(i)
!         end do
!     end function getCovariance

!     function randbetween(a,b) result(random)
!         real, intent(in) :: a, b
!         real :: r, random
!         call RANDOM_NUMBER(r)
!         random = (b-a)*r + a;
!     end

!     function getDistance(a,b) result(distance)
!         real, intent(in) :: a(:), b(:)
!         real :: distance
!         !d = norm((a-b)./fullDist);
!         distance = norm2(a-b);
!     end
    
!     function nearestNeighbour(point, points) result (idx) !O(N)
!         real, intent(in) :: point(:), points(:,:)
!         integer :: idx, i
!         real :: distance, neighbour(size(point)), currentDis
!         distance = huge(distance)
!         idx = -1
!         do i = 1, SIZE(points,2)
!             neighbour = points(:,i);
!             currentDis = getDistance(point, neighbour)
!             if(currentDis < distance) then
!                 distance = currentDis;
!                 idx = i;
!             endif
!         end do
!     end function nearestNeighbour
    
!     subroutine generateGammaMatrix(gammaMatrix, s, x, x1, nSamples) 
!         real, intent(in) :: x(:,:), x1(:,:)
!         integer, intent(in) :: nSamples
!         integer :: numDims, numberOfSamplings, k, d, nearestInX, nearestInX1
!         real :: coordinateRange(SIZE(x,1),2), coord(SIZE(x,1))
!         integer, intent(out) :: gammaMatrix(size(x,2),size(x1,2))
!         integer, intent(out) :: s(size(x1,2))
!         integer :: i
!         numDims = size(x,1)
!         s = 0
!         do d = 1,numDims
!             coordinateRange(d,1) = minval(x1(d,:))
!             coordinateRange(d,2) = maxval(x1(d,:))
!         end do
!         numberOfSamplings = nSamples * size(x,2)
!         gammaMatrix = 0
!         do k = 1, numberOfSamplings
!             coord = 0
!             do d = 1,numDims
!                 coord(d) = randbetween(coordinateRange(d,1),coordinateRange(d,2))
!             end do
!             nearestInX = nearestNeighbour(coord, x)
!             nearestInX1 = nearestNeighbour(coord, x1)
        
!             gammaMatrix(nearestInX,nearestInX1) = gammaMatrix(nearestInX, nearestInX1) + 1
!             s(nearestInX1) = s(nearestInX1) + 1
!         end do

        
        

!         do i = 1,size(s) 
!             if (s(i) == 0) then
                
!                 s(i) = 1
!                 coord = x1(:,i)
!                 nearestInX = nearestNeighbour(coord, x)
!                 gammaMatrix(nearestInX, i) = 1
!             endif
!         end do
        
        
!     end subroutine generateGammaMatrix

 end module linearInterpolation