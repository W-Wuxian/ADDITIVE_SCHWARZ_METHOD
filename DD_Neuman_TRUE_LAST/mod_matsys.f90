!----------------------------------------------------
! CE MODULE CONTIENT LA CONSTRUCTION DE LA MATRICE A
!              ET LE VECTEUR SOURCE F
!----------------------------------------------------

Module mod_matsys
!---modules-------------------
  Use mod_parametres
  Use mod_fonctions
  Use mod_gradconj
  Use mod_constr_mat
!-----------------------------
  Implicit None

Contains
  !--- SUBROUTINE QUI CONSTRUIT LA MATRICE A : UN SEUL APPEL ----------------------

  Subroutine matsys_v2(nnz,Nx_l,Ny_l,AA,IA,JA)
    Integer, Intent(IN) ::nnz,Nx_l,Ny_l
    Real(PR), Dimension(:), Allocatable, Intent(Out) :: AA
    Integer, Dimension(:), Allocatable, Intent(Out) :: IA,JA
    Integer :: k, mul2, d1, d2, mul1
    k = 1; mul2 = Nx_l*Ny_l; mul1 = Ny_l*(Nx_l-1)+1
    ALLOCATE(AA(nnz),IA(nnz),JA(nnz))
    
    !L1-----------------------------------------------------------
    IF(Nx_l > 1)THEN
       if (rang/=0) then
          CALL L1_neuman(Ny_l,nnz,k,AA,IA,JA)
       else
          CALL L1(Ny_l,nnz,k,AA,IA,JA)
       end if
       !CRTL nnz_L1
       IF(k /= crtl_L1_nnz)THEN
          PRINT*,'ERROR L1_nnz','RANG',rang,'k_L1',k,'/=','crtl_L1_nnz',crtl_L1_nnz
          AA=0._PR;IA=0;JA=0
       END IF
    ELSE IF(Nx_l == 1)THEN
       CALL L1_NXL_eqv_1(Ny_l,nnz,k,AA,IA,JA)
       IF(k /= crtl_L1_nnz_Nxl_EQV_1)THEN
          PRINT*,'ERROR L1_nnz','RANG',rang,'k_L1',k,'/=','crtl_L1_nnz_Nxl_EQV_1',crtl_L1_nnz_Nxl_EQV_1
          AA=0._PR;IA=0;JA=0
       END IF
    END IF
    
    !L2-----------------------------------------------------------
    IF(Nx_l>2)THEN
       d1 = Ny_l+1; d2 = 2*Ny_l
       CALL  L2(Nx_l,Ny_l,nnz,d1,d2,k,AA,IA,JA)
       IF(k /= crtl_L1_nnz + crtl_L2_nnz)THEN
          PRINT*,'ERROR L2_nnz','RANG',rang,'k_L2',k-crtl_L1_nnz,'/=','crtl_L2_nnz',crtl_L2_nnz
          AA=0._PR;IA=0;JA=0
       END IF
    END IF
    !L3-----------------------------------------------------------
    IF(Nx_l>1)THEN
       if (rang/=Np-1) then
          CALL L3_neuman(mul1,mul2,Ny_l,nnz,k,AA,IA,JA)
       else
          CALL L3(mul1,mul2,Ny_l,nnz,k,AA,IA,JA)
       end if
       PRINT*,'rang',rang,'k',k,'FIN L3'
       IF(k /= sum_crtl_L_nnz)THEN
          PRINT*,'ERROR L3_nnz','RANG',rang,'k_L3',k-crtl_L1_nnz-crtl_L2_nnz,'/=','crtl_L3_nnz',crtl_L3_nnz
          AA=0._PR;IA=0;JA=0
       END IF
    END IF
    
    
  END Subroutine matsys_v2


!--------SUBROUTINE DU VECTEUR SOURCE EN FONCTION DE L'ITERATION N ET DU VECTEUR U--------
Subroutine vectsource(CT,U,S1,S2,X,Y,T,F)
  Implicit None
  !---variables----------------------------------
  Integer, Intent(In) :: CT, S1, S2
  Real(PR), Dimension(it1:itN), Intent(In) :: U
  Real(PR), Dimension(LBX:UBX), Intent(IN) :: X
  Real(PR), Dimension(0:Ny_g+1), Intent(IN) :: Y !BECAUSE DD 1D => Ny_l == Ny_g
  Real(PR), Intent(IN) :: T
  Real(PR), Dimension(it1:itN), Intent(INOUT) :: F
  Integer :: i,j,k
  
  DO i = S1, S2
     DO j = 1, Ny_g
        k = j + (i-1)*Ny_g
        F(k) = dt * f1(CT,X(i),Y(j),T)
        IF(i == 1)THEN
           F(k) = F(k) - h1(CT,X(i-1),Y(j),T)*beta
        ELSE IF(i == Nx_g)THEN
           F(k) = F(k) - h1(CT,X(i+1),Y(j),T)*beta
        END IF
        IF(j == 1)THEN
           F(k) = F(k) - g1(CT,X(i),Y(j-1),T)*gamma
        ELSE IF(j == Ny_g)THEN
           F(k) = F(k) - g1(CT,X(i),Y(j+1),T)*gamma
        END IF
     END DO
  END DO
  F = F +U
End Subroutine vectsource



Subroutine vectsource_FULL(CT,U,UG,UD,S1,S2,X,Y,T,F)
  Implicit None
  !---variables----------------------------------
  Integer, Intent(In) :: CT, S1, S2
  Real(PR), Dimension(it1:itN), Intent(In) :: U
  Real(PR), Dimension(LBUG:UBUG), Intent(In) :: UG
  Real(PR), Dimension(LBUD:UBUD), Intent(In) :: UD
  Real(PR), Dimension(LBX:UBX), Intent(IN) :: X
  Real(PR), Dimension(0:Ny_g+1), Intent(IN) :: Y !BECAUSE DD 1D => Ny_l == Ny_g
  Real(PR), Intent(IN) :: T
  Real(PR), Dimension(it1:itN), Intent(INOUT) :: F
  Integer :: i,j,k
  
  DO i = S1, S2
     DO j = 1, Ny_g
        k = j + (i-1)*Ny_g
        F(k) = dt * f1(CT,X(i),Y(j),T)
        IF(i == 1)THEN
           F(k) = F(k) - h1(CT,X(i-1),Y(j),T)*beta
        ELSE IF(i == S1 )THEN !.AND. rang /= 0)THEN not needed with the current order if... else if
           F(k) = F(k) + (D*dt/dx)*UG(k)*c_neuman + beta*( UG(k+Ny_g)-UG(k) )  
        ELSE IF(i == Nx_g)THEN
           F(k) = F(k) - h1(CT,X(i+1),Y(j),T)*beta
        ELSE IF(i == S2)THEN !.AND. rang /= Np-1)THEN not needed with the current order else if... else if
           F(k) = F(k) + (D*dt/dx)*UD(k)*c_neuman + beta*( UD(k-Ny_g)-UD(k) )
        END IF
        IF(j == 1)THEN
           F(k) = F(k) - g1(CT,X(i),Y(j-1),T)*gamma
        ELSE IF(j == Ny_g)THEN
           F(k) = F(k) - g1(CT,X(i),Y(j+1),T)*gamma
        END IF
     END DO
  END DO
  F = F +U
End Subroutine vectsource_FULL

End Module mod_matsys






!!$DO d = 1, Ny_l
!!$       IF(d == 1)THEN
!!$          AA(k)=alpha; IA(k)=d; JA(k)=d
!!$          DO hit = d+1, d+Ny_l, Ny_l-1
!!$             IF(hit == d+1)THEN
!!$                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=hit
!!$             ELSE
!!$                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=hit
!!$             END IF
!!$          END DO
!!$       ELSE IF(d>1 .AND. d<Ny_l)THEN
!!$          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
!!$          k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
!!$          k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
!!$          k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
!!$       ELSE IF(d == Ny_l)THEN
!!$           DO hit = d-1, d
!!$             IF(hit == d-1)THEN
!!$                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=hit
!!$             ELSE
!!$                k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=hit
!!$             END IF
!!$             k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
!!$          END DO
!!$       END IF
!!$    END DO







!!$IF(Nx_l>2)THEN
!!$       d1 = Ny_l+1; d2 = 2*Ny_l
!!$       DO i = 1, Nx_l-2
!!$          DO d = d1,d2
!!$             IF(d == d1)THEN
!!$                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
!!$                k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
!!$                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
!!$                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
!!$             ELSE IF(d>d1 .AND. d<d2)THEN
!!$                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
!!$                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
!!$                k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
!!$                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d+1
!!$                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
!!$             ELSE IF(d == d2)THEN
!!$                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d-Ny_l
!!$                k=k+1; AA(k)=gamma; IA(k)=d; JA(k)=d-1
!!$                k=k+1; AA(k)=alpha; IA(k)=d; JA(k)=d
!!$                k=k+1; AA(k)=beta; IA(k)=d; JA(k)=d+Ny_l
!!$             END IF
!!$          END DO
!!$          d1 = d2+1; d2=d2+Ny_l
!!$       END DO
!!$    END IF


























!!$Subroutine matsys(nnz,AA,IA,JA)
!!$  Integer, Intent(In) :: nnz
!!$  Real(PR), Dimension(:), Allocatable, Intent(Out) :: AA
!!$  Integer, Dimension(:), Allocatable, Intent(Out) :: IA,JA
!!$  Integer :: i,j,k
!!$  k=1
!!$  Allocate(AA(nnz),IA(nnz),JA(nnz))
!!$  Do i=1,na
!!$    Do j=1,na
!!$      If (i==j) Then !diagonale principale
!!$        AA(k)=alpha
!!$        IA(k)=i
!!$        JA(k)=j
!!$        k=k+1
!!$      ElseIf ((i==j-1) .and. (modulo(i,nx)/=0)) Then !diagonale sup
!!$        AA(k)=beta
!!$        IA(k)=i
!!$        JA(k)=j
!!$        k=k+1
!!$      ElseIf ((i==j+1) .and. (i/=1) .and. (modulo(j,ny)/=0)) Then !diagonale inf
!!$        AA(k)=beta
!!$        IA(k)=i
!!$        JA(k)=j
!!$        k=k+1
!!$      ElseIf (j==ny+i) Then !diagonale la plus a droite
!!$        AA(k)=gamma
!!$        IA(k)=i
!!$        JA(k)=j
!!$        k=k+1
!!$      ElseIf (i==j+nx) Then !diagonale la plus a gauche
!!$        AA(k)=gamma
!!$        IA(k)=i
!!$        JA(k)=j
!!$        k=k+1
!!$      End If
!!$    End Do
!!$  End Do
!!$End Subroutine matsys
