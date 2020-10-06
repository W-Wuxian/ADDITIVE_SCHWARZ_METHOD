module ModF
  use mpi

  implicit none
INTEGER,parameter::RK_F=8
! ModF CONTIENT:
!SUBROUTINES PERMETTANT DE CALCULER LES FONCTIONS SOURCES,
!LES FONCTIONS DE BORDS 
contains


  SUBROUTINE F1(it1,itN,S1,S2,INTER_Y,GAP,X,Y,F)
    INTEGER,INTENT(IN)::it1,itN,S1,S2,INTER_Y
    INTEGER,DIMENSION(S1:S2,1:2),INTENT(IN)::GAP
    REAL(RK_F),DIMENSION(S1:S2),INTENT(IN)::X
    REAL(RK_F),DIMENSION(1:INTER_Y),INTENT(IN)::Y
    REAL(RK_F),DIMENSION(it1:itN),INTENT(OUT)::F
    INTEGER::i,j,k
    REAL(RK_F)::FX
    DO i=S1,S2!CHARGE LOCALE EN X
       FX=X(i)-X(i)**2!POUR i FIXER, ON PEUT NE CALCULER FX QU'UNE SEULE FOIS
       DO j=GAP(i,1),GAP(i,2)!CHARGE LOCALE EN Y
          k=j+(i-1)*INTER_Y!CHARGE GLOBALE
          F(k)=2*(Y(j)-Y(j)**2)+2*FX!CALCUL FONCTION SOURCE
       END DO
    END DO
  END SUBROUTINE F1


  SUBROUTINE F2(it1,itN,S1,S2,INTER_Y,GAP,X,Y,F)
    INTEGER,INTENT(IN)::it1,itN,S1,S2,INTER_Y
    INTEGER,DIMENSION(S1:S2,1:2),INTENT(IN)::GAP
    REAL(RK_F),DIMENSION(S1:S2),INTENT(IN)::X
    REAL(RK_F),DIMENSION(1:INTER_Y),INTENT(IN)::Y
    REAL(RK_F),DIMENSION(it1:itN),INTENT(OUT)::F
    INTEGER::i,j,k
    REAL(RK_F)::FX
    DO i=S1,S2
       FX=sin(X(i))!POUR i FIXER, ON PEUT NE CALCULER FX QU'UNE SEULE FOIS
       DO j=GAP(i,1),GAP(i,2)
          k=j+(i-1)*INTER_Y
          F(k)=cos(Y(j))+FX
       END DO
    END DO
  END SUBROUTINE F2
  
  SUBROUTINE F3_FIXE(Lx,Ly,it1,itN,S1,S2,INTER_Y,GAP,X,Y,F)!ON NE STOCKE DANS LE VECTEUR F,
    !QUE LA PARTIE FIXE DE F3
    REAL(RK_F),INTENT(IN)::Lx,Ly
    INTEGER,INTENT(IN)::it1,itN,S1,S2,INTER_Y
    INTEGER,DIMENSION(S1:S2,1:2),INTENT(IN)::GAP
    REAL(RK_F),DIMENSION(S1:S2),INTENT(IN)::X
    REAL(RK_F),DIMENSION(1:INTER_Y),INTENT(IN)::Y
    REAL(RK_F),DIMENSION(it1:itN),INTENT(OUT)::F
    INTEGER::i,j,k
    REAL(RK_F)::FX
    DO i=S1,S2
       FX=exp(-(X(i)-0.50d0*Lx)**2)!POUR i FIXER, ON PEUT NE CALCULER FX QU'UNE SEULE FOIS
       DO j=GAP(i,1),GAP(i,2)
          k=j+(i-1)*INTER_Y
          F(k)=FX*exp(-(Y(j)-0.50d0*Ly)**2)
       END DO
    END DO
  END SUBROUTINE F3_FIXE


  SUBROUTINE WR_F(rang,it1,itN,INTER_Y,F,S1,S2,GAP,X,Y)! ECRITUE DE F DANS LE REPERTOIRE: F  
    INTEGER,INTENT(IN)::rang,it1,itN,S1,S2,INTER_Y
    INTEGER,DIMENSION(S1:S2,1:2),INTENT(IN)::GAP
    REAL(RK_F),DIMENSION(S1:S2),INTENT(IN)::X
    REAL(RK_F),DIMENSION(1:INTER_Y),INTENT(IN)::Y
    REAL(RK_F),DIMENSION(it1:itN),INTENT(IN)::F
    INTEGER::i,j,k
    CHARACTER(LEN=3)::rank
    WRITE(rank,fmt='(1I3)')rang
    OPEN(rang+10,file='F/F_'//trim(adjustl(rank))//'.dat')
    DO i=S1,S2
       DO j=GAP(i,1),GAP(i,2)
          k=j+(i-1)*INTER_Y
          WRITE(10+rang,*)X(i),Y(j),F(k)
       END DO
    END DO
    CLOSE(10+rang)
  END SUBROUTINE WR_F

  SUBROUTINE BC(AP,AL,userchoice,S1,S2,GAP,it1,itN,INTER_Y,INTER_X&
       ,X,Y,DX,DY,Lx,Ly,VECT)!ON CALCUL LES CONDITIONS LIMITES(BC)
    REAL(RK_F),INTENT(IN)::AP,AL!LES COEFF DE LA MATRICE
    INTEGER,INTENT(IN)::userchoice,S1,S2,it1,itN,INTER_X,INTER_Y
    INTEGER,DIMENSION(S1:S2,1:2),INTENT(IN)::GAP
    REAL(RK_F),INTENT(IN)::DX,DY,Lx,Ly
    REAL(RK_F),DIMENSION(S1:S2),INTENT(IN)::X
    REAL(RK_F),DIMENSION(1:INTER_Y),INTENT(IN)::Y
    REAL(RK_F),DIMENSION(it1:itN),INTENT(INOUT)::VECT
    integer::i,j,k
    real(RK_F)::G,H
    SELECT CASE(userchoice)
    CASE(2)!F2
       DO i=S1,S2
          DO j=GAP(i,1),GAP(i,2)
              k=j+(i-1)*INTER_Y
              IF(i==1)THEN!BC H
              VECT(k)=VECT(k)+AL*cos(Y(j))
             ELSE IF(i==INTER_X)THEN
                VECT(k)=VECT(k)+AL*(cos(Y(j))+sin(Lx))
             END IF
             IF(j==1)THEN!BC G
                VECT(k)=VECT(k)+AP*(1.0d0+sin(X(i)))
             ELSE IF(j==INTER_Y)THEN
                VECT(k)=VECT(k)+AP*(cos(Ly)+sin(X(i)))
             END IF
          END DO
       END DO
    CASE(3)!F3
       IF(S1==1)THEN!BC
          i=S1
          DO j=GAP(i,1),GAP(i,2)
             k=j+(i-1)*INTER_Y;VECT(k)=VECT(k)+AL !CAR:G=0 ET H=1
          END DO
       END IF
       IF(S2==INTER_X)THEN!BC
          i=S2
          DO j=GAP(i,1),GAP(i,2)
             k=j+(i-1)*INTER_Y;VECT(k)=VECT(k)+AL !CAR:G=0 ET H=1
          END DO
       END IF
    END SELECT
  END SUBROUTINE BC
  
end module ModF
