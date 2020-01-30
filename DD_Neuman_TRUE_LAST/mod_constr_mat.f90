module mod_constr_mat
  use mod_parametres

  implicit none

!-->set a_neuman=1.0_PR et b_neuman=0.0000001_PR pour retrouver les mêmes ordres de grandeurs des erreurs L2 qu'avec le sequentiel ou dirichlet 
Real(PR), Parameter :: a_neuman = 0.5_PR  
Real(PR), Parameter :: b_neuman = 0.5_PR
REAL(PR), Parameter :: c_neuman = a_neuman/b_neuman 
Real(PR), Parameter :: alpha_neuman_H = 1._PR  -2._PR*beta -2._PR*gamma +((D*dt*a_neuman)/(dx*b_neuman)) +beta
!  Les coefficients alpha changent pour les Ny première
!  et/ou dernières lignes des matrices dans le cas des CLs de Neuman
Real(PR), Parameter :: alpha_neuman_B = 1._PR -2._PR*beta -2._PR*gamma +((D*dt*a_neuman)/(dx*b_neuman)) +beta

contains

  Subroutine L1_NXL_eqv_1(Ny_l,nnz,k,AA,IA,JA)
    INTEGER,INTENT(IN)::Ny_l,nnz
    INTEGER,INTENT(INOUT)::k
    Real(PR), Dimension(nnz), Intent(Out) :: AA
    Integer, Dimension(nnz), Intent(Out) :: IA,JA
    Integer::d
    !ONLY WHEN A PROCES  GOT Nx_l==1
    DO d = 1, Ny_l
       IF(d == 1)THEN
          AA(k)=alpha; IA(k)=d; JA(k)=d
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
       ELSE IF(d>1 .AND. d<Ny_l)THEN
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
          k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
       ELSE IF(d == Ny_l)THEN
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
          k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
       END IF
    END DO
  End Subroutine L1_NXL_eqv_1
  
  Subroutine L1(Ny_l,nnz,k,AA,IA,JA)
    INTEGER,INTENT(IN)::Ny_l,nnz
    INTEGER,INTENT(INOUT)::k
    Real(PR), Dimension(nnz), Intent(Out) :: AA
    Integer, Dimension(nnz), Intent(Out) :: IA,JA
    Integer::d,hit
    DO d = 1, Ny_l
       IF(d == 1)THEN
          AA(k)=alpha; IA(k)=d; JA(k)=d
          DO hit = d+1, d+Ny_l, Ny_l-1
             IF(hit == d+1)THEN
                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=hit
             ELSE
                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=hit
             END IF
          END DO
       ELSE IF(d>1 .AND. d<Ny_l)THEN
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
          k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
          k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
       ELSE IF(d == Ny_l)THEN
          DO hit = d-1, d
             IF(hit == d-1)THEN
                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=hit
             ELSE 
                k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=hit
             END IF
          END DO
          k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
       END IF
    END DO
  End Subroutine L1

  SUBROUTINE TRACK(k,d)
    INTEGER,INTENT(INOUT)::k,d
    IF(rang==0)THEN
       PRINT*,'RANG',rang,'TRACK k,d',k,d
    END IF
  END SUBROUTINE TRACK


  Subroutine L2(Nx_l,Ny_l,nnz,d1,d2,k,AA,IA,JA)
    Integer, Intent(IN)    :: Nx_l, Ny_l, nnz
    Integer, Intent(INOUT) :: d1, d2, k
    Real(PR), Dimension(nnz), Intent(Out) :: AA
    Integer, Dimension(nnz), Intent(Out) :: IA,JA
    INTEGER :: i, d
    DO i = 1, Nx_l-2
       DO d = d1,d2
          IF(d == d1)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
          ELSE IF(d>d1 .AND. d<d2)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
             k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
          ELSE IF(d == d2)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
             k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
          END IF
       END DO
       d1 = d2+1; d2=d2+Ny_l
    END DO
  End Subroutine L2


  Subroutine L3(mul1,mul2,Ny_l,nnz,k,AA,IA,JA)
    Integer, Intent(IN) :: mul1, mul2, Ny_l, nnz
    Integer, Intent(INOUT) :: k
    Real(PR), Dimension(nnz), Intent(Out) :: AA
    Integer, Dimension(nnz), Intent(Out) :: IA,JA
    INTEGER :: d
    DO d = mul1, mul2
          IF(d == mul1)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
          ELSE IF(d>mul1 .AND. d<mul2)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
             k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
          ELSE IF(d == mul2)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
             k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
          END IF
       END DO
  End Subroutine L3

  Subroutine L1_Neuman(Ny_l,nnz,k,AA,IA,JA)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::Ny_l,nnz
    INTEGER,INTENT(INOUT)::k
    Real(PR), Dimension(nnz), Intent(Out) :: AA
    Integer, Dimension(nnz), Intent(Out) :: IA,JA
    Integer::d,hit
    
    DO d = 1, Ny_l
       IF(d == 1)THEN
          AA(k)=alpha_neuman_H
          IA(k)=d
          JA(k)=d
          DO hit = d+1, d+Ny_l, Ny_l-1
             IF(hit == d+1)THEN
                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=hit
             ELSE
                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=hit
             END IF
          END DO
       ELSE IF(d>1 .AND. d<Ny_l)THEN
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
          k=k+1; AA(k)=alpha_neuman_H; IA(k)=d; JA(k)=d
          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
          k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
       ELSE IF(d == Ny_l)THEN
          DO hit = d-1, d
             IF(hit == d-1)THEN
                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=hit
             ELSE 
                k=k+1; AA(k)=alpha_neuman_H; IA(k)=d; JA(k)=hit
             END IF
          END DO
          k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
       END IF
    END DO
  End Subroutine L1_Neuman

  Subroutine L3_neuman(mul1,mul2,Ny_l,nnz,k,AA,IA,JA)
    IMPLICIT NONE

    Integer, Intent(IN) :: mul1, mul2, Ny_l, nnz
    Integer, Intent(INOUT) :: k
    Real(PR), Dimension(nnz), Intent(Out) :: AA
    Integer, Dimension(nnz), Intent(Out) :: IA,JA
    INTEGER :: d
    DO d = mul1, mul2
          IF(d == mul1)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=alpha_neuman_B; IA(k)=d; JA(k)=d
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
          ELSE IF(d>mul1 .AND. d<mul2)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
             k=k+1; AA(k)=alpha_neuman_B; IA(k)=d; JA(k)=d
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
          ELSE IF(d == mul2)THEN
             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
             k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
             k=k+1; AA(k)=alpha_neuman_B; IA(k)=d; JA(k)=d
          END IF
       END DO
  End Subroutine L3_neuman

end module mod_constr_mat
