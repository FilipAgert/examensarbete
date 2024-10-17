MODULE potsurf
  
   ! --------------------------------------------------------------------
   ! Data types
   ! --------------------------------------------------------------------
   
   INTEGER(4), PARAMETER :: r_kind=8 ! Number of bytes for real FP numbers
   INTEGER(4), PARAMETER :: i_kind=4 ! Number of bytes for integer numbers
   
   ! --------------------------------------------------------------------
   ! Potential surface
   ! --------------------------------------------------------------------
   
   ! Limits of the parameter grid
   ! The meaning of the indexes are on page 8 in PRC84, 034613 (2011)
   ! Details about the grid can be found in PRC79, 064304 (2009)
   INTEGER(kind=i_kind), parameter :: MIN_II = 1, MAX_II = 53 !Elongation Q
   INTEGER(kind=i_kind), parameter :: MIN_JJ = 1, MAX_JJ = 15 !Neck
   INTEGER(kind=i_kind), parameter :: MIN_KK = 1, MAX_KK = 15 !Left fragment deformation
   INTEGER(kind=i_kind), parameter :: MIN_LL = 1, MAX_LL = 15 !Right fragment deformation
   INTEGER(kind=i_kind), parameter :: MIN_MM = -30, MAX_MM = 30 !Mass/reflection assymetry
   
   REAL(kind=r_kind) :: Epot(MIN_II:MAX_II,MIN_JJ:MAX_JJ,MIN_KK:MAX_KK,MIN_LL:MAX_LL,MIN_MM:MAX_MM)         ! Potential energy
   REAL(kind=r_kind) :: Emac(MIN_II:MAX_II,MIN_JJ:MAX_JJ,MIN_KK:MAX_KK,MIN_LL:MAX_LL,MIN_MM:MAX_MM)         ! Macroscopic energy
   REAL(kind=r_kind) :: Rneck(MIN_II:MAX_II,MIN_JJ:MAX_JJ,MIN_KK:MAX_KK,MIN_LL:MAX_LL,MIN_MM:MAX_MM)         ! Neck radii
 
   
   
 CONTAINS
 
 
   REAL(kind=r_kind) FUNCTION get_q2_value(AZ,AA,I_val) RESULT(q2_val)
 
     ! Returns the reduced quadrupole moment q2
 
     ! Input
     ! AZ: Proton number
     ! AA: Mass number
     ! I_val: Grid index I
     
     IMPLICIT NONE
     INTEGER(kind=i_kind), intent(in) :: AZ,AA
     INTEGER(kind=i_kind), intent(in) :: I_val
     REAL(kind=r_kind), parameter :: Pi = ACOS(-1.0_r_kind)
     REAL(kind=r_kind) :: q2_scale
     REAL(kind=r_kind) :: Q2factor !scaling factor for Q2values
     INTEGER(kind=i_kind), parameter :: Aref = 240
     INTEGER(kind=i_kind), parameter :: Zref = 94
     REAL(kind=r_kind), parameter :: onethird = 1.0_r_kind/3.0_r_kind
     REAL(kind=r_kind), parameter :: twothird = 2.0_r_kind/3.0_r_kind
     
     ! Quadrupole deformation values for 240Pu - corresponding to index II.
     REAL(kind=r_kind), parameter :: Q2vals(53) = (/4.81253,   6.82805,   8.90879,  11.06605,  13.31188,  15.65945,  18.12340,&
       20.72011,  23.46814,  26.38871,  29.50621,  32.84902,  36.45030,  38.39975,  40.34920,  42.47073,  44.59226,  46.91380,&
         49.23533,  51.79066,  54.34598,  57.17638,  60.00677,  63.12872,  66.25067,  69.48992,  72.72917,  76.04370,  79.35822,&
           86.12054,  92.99791,  99.97944, 107.05082, 114.19157, 121.42655, 128.69212, 136.03868, 143.41393, 150.81228, 158.20718,&
            165.78249, 173.09174, 180.81321, 188.18743, 195.65203, 203.06961, 210.66009, 218.12864, 225.62603, 233.36657,&
             240.74626, 247.79498, 254.57527/)
 
     Q2factor = REAL(AZ,kind=r_kind)/REAL(Zref,kind=r_kind)*( REAL(AA,kind=r_kind)/REAL(Aref,kind=r_kind) )**twothird
     
     q2_scale = 0.01_r_kind*3*AZ*1.2_r_kind**2*REAL(AA,kind=r_kind)**(2.0_r_kind/3.0_r_kind)/(4.0_r_kind*Pi) ! 0.01 factor as Q2 is in barns=100fm^2
 
     q2_val = Q2factor*Q2vals(I_val)/q2_scale
     
   END FUNCTION get_q2_value
 
 
   
   
   
   SUBROUTINE read_5D_epot(AZ,AA,filename_pot5D)
     IMPLICIT NONE
     
     !Set up the matrix of energies,neck radii using the read values. 
     !The matrix is set up so that points obtained by rotating 
     !the nucleus by pi around the axis perpendiucalr to the axial symmetry axis are included.
     
     INTEGER(kind=i_kind), intent(in) :: AZ,AA
     CHARACTER(*), intent(in) :: filename_pot5D
     REAL(kind=r_kind) :: Evals(MIN_II:MAX_II,MIN_JJ:MAX_JJ,MIN_KK:MAX_KK,MIN_LL:MAX_LL,1:MAX_MM+2)
     INTEGER(kind=i_kind) :: iounit
     INTEGER(kind=i_kind) :: II, JJ, KK, LL, MM
     REAL(kind=r_kind) :: Eval
     INTEGER(kind=i_kind) :: Z_val, A_val
     INTEGER(kind=i_kind) :: as, eL, eR
     INTEGER :: point
     
     !set all points to unphysical values, should be overwritten when file is read.
     Evals = 1000
     
     WRITE(*,*)''
     WRITE(*,*) 'Reading energies from file ',filename_pot5D
     OPEN(file=filename_pot5D, status='old', action='read',newunit=iounit)
     READ(iounit,*) Z_val, A_val
     
     IF(Z_val /= AZ) THEN
        WRITE(*,*) 'WARNING: Z in inputfile', filename_pot5D, 'is Z= ',Z_val,'does not match requested AZ= ',AZ
        !STOP
     END IF
     IF(A_val /= AA) THEN
        WRITE(*,*) 'WARNING: A in inputfile', filename_pot5D, 'is A= ',A_val,'does not match requested AA= ',AA
        !STOP
     END IF
     
     point = 0
     
     DO II = MIN_II,MAX_II
        DO JJ = MIN_JJ,MAX_JJ
           DO KK = MIN_KK,MAX_KK
              DO LL = MIN_LL,MAX_LL
                 DO MM = 1,MAX_MM+2 !due to symmetry, PM has one extra value at MM = 1 which is not used
                    
                    point = point + 1
                    
                    READ(iounit,'(f9.3)') Eval
                    Evals(II,JJ,KK,LL,MM) = Eval
                    
                 END DO
              END DO
           END DO
        END DO
     END DO
     
     CLOSE(iounit)
     
 
     ! Shift so that MM=0 is symmetric mass
     Epot(:,:,:,:,0:MAX_MM) = Evals(:,:,:,:,2:MAX_MM+2)
     
     
   DO as = 1,MAX_MM
      DO eL = MIN_KK,MAX_KK
         DO eR = MIN_LL,MAX_LL
            Epot(:,:,eR,eL,-as) = Epot(:,:,eL,eR,as)
         END DO
      END DO
   END DO
     
     WRITE(*,*) 'Done reading, read ',point,' points'
     
   END SUBROUTINE read_5D_epot
 
 
 
 
 
 
   SUBROUTINE read_5d_emac(AZ,AA,filename_emac5D)
     IMPLICIT NONE
     
     INTEGER(kind=i_kind), intent(in) :: AZ,AA
     CHARACTER(*), intent(in) :: filename_emac5D 
     INTEGER(kind=i_kind) :: iounit
     INTEGER(kind=i_kind) :: II, JJ, KK, LL, MM
     INTEGER(kind=i_kind) :: as, eL, eR
     INTEGER(kind=i_kind) :: Z_val, A_val
     REAL(kind=r_kind) :: Eval
     REAL(kind=r_kind) :: Emacvals(MIN_II:MAX_II,MIN_JJ:MAX_JJ,MIN_KK:MAX_KK,MIN_LL:MAX_LL,1:MAX_MM+2)
     INTEGER :: point
     
     WRITE(*,*)''
     WRITE(*,*) 'Reading macroscopic energies from file ',filename_emac5D
     OPEN(file=filename_emac5D, status='old', action='read',newunit=iounit)
     
     READ(iounit,*) Z_val, A_val
     
     IF(Z_val /= AZ) THEN
        WRITE(*,*) 'WARNING: Z in inputfile', filename_emac5D, 'is Z= ',Z_val,'does not match requested AZ= ',AZ
        !STOP
     END IF
     IF(A_val /= AA) THEN
        WRITE(*,*) 'WARNING: A in inputfile', filename_emac5D, 'is A= ',A_val,'does not match requested AA= ',AA
        !STOP
     END IF
     
     point = 0
     
     DO II = MIN_II,MAX_II
        DO JJ = MIN_JJ,MAX_JJ
           DO KK = MIN_KK,MAX_KK
              DO LL = MIN_LL,MAX_LL
                 DO MM = 1,MAX_MM+2 !due to symmetry, PM has one extra value at MM = 1 which is not used
                    
                    point = point + 1
                    
                    READ(iounit,'(f9.3)') Eval
                    Emacvals(II,JJ,KK,LL,MM) = Eval
                    
                 END DO
              END DO
           END DO
        END DO
     END DO
     
     CLOSE(iounit)
     
     Emac(:,:,:,:,0:MAX_MM) = Emacvals(:,:,:,:,2:MAX_MM+2)
     
   DO as = 1,MAX_MM
      DO eL = MIN_KK,MAX_KK
         DO eR = MIN_LL,MAX_LL
            Emac(:,:,eR,eL,-as) = Emac(:,:,eL,eR,as)           
         END DO
      END DO
   END DO
     
     WRITE(*,*) 'Done reading, read ',point,' points'
     
   END SUBROUTINE read_5d_emac
   
 
 
   SUBROUTINE read_5d_rneck(AZ,AA,filename_rneck5D)
     IMPLICIT NONE
     
     INTEGER(kind=i_kind), intent(in) :: AZ,AA
     CHARACTER(*), intent(in) :: filename_rneck5D 
     INTEGER(kind=i_kind) :: iounit
     INTEGER(kind=i_kind) :: II, JJ, KK, LL, MM
     INTEGER(kind=i_kind) :: as, eL, eR
     INTEGER(kind=i_kind) :: Z_val, A_val
     REAL(kind=r_kind) :: rneckval
     REAL(kind=r_kind) :: Rneckvals(MIN_II:MAX_II,MIN_JJ:MAX_JJ,MIN_KK:MAX_KK,MIN_LL:MAX_LL,1:MAX_MM+2)
     INTEGER :: point
     REAL(kind=r_kind) :: rscaling
     REAL(kind=r_kind) :: r0 = 1.2_r_kind
 
     Rneckvals = 0.0
     
     WRITE(*,*)''
     WRITE(*,*) 'Reading neck radii from file ',filename_rneck5D
     OPEN(file=filename_rneck5D, status='old', action='read',newunit=iounit)
     
     point = 0
     
     DO II = MIN_II,MAX_II
        DO JJ = MIN_JJ,MAX_JJ
           DO KK = MIN_KK,MAX_KK
              DO LL = MIN_LL,MAX_LL
                 DO MM = 1,MAX_MM+2 !due to symmetry, PM has one extra value at MM = 1 which is not used
 
                    point = point + 1
                    
                    READ(iounit,*) rneckval
                    Rneckvals(II,JJ,KK,LL,MM) = rneckval
 
                 END DO
              END DO
           END DO
        END DO
     END DO
     
     CLOSE(iounit)
 
     !rescales the neck radius for the current nucleus
     rscaling = r0*REAL(AA,kind=r_kind)**(1.0_r_kind/3.0_r_kind)
     
     Rneck(:,:,:,:,0:MAX_MM) = rscaling*Rneckvals(:,:,:,:,2:MAX_MM+2)
     
     DO as = 1,MAX_MM
        DO eL = MIN_KK,MAX_KK
           DO eR = MIN_LL,MAX_LL
              Rneck(:,:,eR,eL,-as) = Rneck(:,:,eL,eR,as)
           END DO
        END DO
     END DO
     
     WRITE(*,*) 'Done reading, read ',point,' points'
     
   END SUBROUTINE read_5d_rneck
 
 
 
   REAL(kind=r_kind) FUNCTION transition_probability(AZ,AA,Etot,I1,J1,K1,L1,M1,I2,J2,K2,L2,M2) RESULT(p)
     ! recipe in PRC 88, 064606 (2013)
 
     ! Calculates the transition probability
     ! between two shape 1 and shape 2
     
     IMPLICIT NONE
 
     INTEGER(kind=i_kind), intent(in) :: AZ, AA
     REAL(kind=r_kind), intent(in) :: Etot
     INTEGER(kind=i_kind), intent(in) :: I1,J1,K1,L1,M1
     INTEGER(kind=i_kind), intent(in) :: I2,J2,K2,L2,M2
     REAL(kind=r_kind) :: aparam
     REAL(kind=r_kind) :: Eexc_1, Eexc_2
     REAL(kind=r_kind) :: Emac_1, Emac_2
     REAL(kind=r_kind) :: Epot_1, Epot_2
     REAL(kind=r_kind) :: Ush_1, Ush_2
     REAL(kind=r_kind) :: Ueff_1, Ueff_2
     REAL(kind=r_kind) :: Eeff_1, Teff_1      
     REAL(kind=r_kind) :: fg_e0param = 8.0_r_kind     ! e0 in expression for Fermi-gas level density parameter a = A/e0
 
     ! Level-density parameter
     aparam = REAL(AA,r_kind)/fg_e0param
 
     ! Potential energy
     Epot_1 = Epot(I1,J1,K1,L1,M1)
     Epot_2 = Epot(I2,J2,K2,L2,M2)
 
     ! Macroscopic energy
     Emac_1 = Emac(I1,J1,K1,L1,M1)
     Emac_2 = Emac(I2,J2,K2,L2,M2)
 
     ! Excitation energy
     Eexc_1 = Etot - Epot_1
     Eexc_2 = Etot - Epot_2
 
     IF(Eexc_2 <= 1.0e-15_r_kind) THEN !no tunneling
        p = 0.0_r_kind
        RETURN
     END IF
 
     ! Shell correction energy
     Ush_1 = Epot_1 - Emac_1 
     Ush_2 = Epot_2 - Emac_2 
 
     !Eq. (4) in PRC88, 064606
     Ueff_1 = Emac_1 + S_function(Eexc_1)*Ush_1
     Ueff_2 = Emac_2 + S_function(Eexc_2)*Ush_2
 
     !Eq. (5)
     Eeff_1 = Eexc_1 + (1.0_r_kind - S_function(Eexc_1))*Ush_1
 
     !Eq. (15) 
     Teff_1 = sqrt ( Eeff_1 / aparam )    
     
     !Eq. (18)
     p =  MIN( EXP( -(Ueff_2-Ueff_1)/Teff_1), 1.0_r_kind ) 
 
   END FUNCTION transition_probability
 
 
 
 
   
   REAL(kind=r_kind) FUNCTION S_function(Eexc)
     IMPLICIT NONE
 
     REAL(kind=r_kind), intent(in) :: Eexc
     REAL(kind=r_kind), parameter :: E0 = 15.0
     REAL(kind=r_kind), parameter :: E1 = 20.0
     REAL(kind=r_kind), parameter :: c = EXP(E1/E0)
  
     !Eq. (19) in PRC88, 064606
     S_function = (1.0_r_kind + c)/(Exp(Eexc/E0) + c)
 
   END FUNCTION S_function   
 
 
 
END MODULE potsurf