module ModSOL
  use mpi

  IMPLICIT NONE
  INTEGER,PARAMETER::RK_SOL=8
!*** MODULE ModSOL POUR DONNER LA SOLUTION EXACTE ET
!*** POUR GERER L'ENREGISTREMENT DE LA SOLUTION DANS DES FICHIERS
!*** SUBROUTINE POUR APPELLER UN SCRIPT GNUPLOT POUR LA VISUALISATION


CONTAINS

  SUBROUTINE UEXACT(U,userchoice,it1,itN,INTER_Y,INTER_X,X,Y,S1,S2,GAP,stateinfo,rang)
    INTEGER,INTENT(IN)::userchoice,it1,itN,INTER_Y,INTER_X,S1,S2,rang
    INTEGER,DIMENSION(S1:S2,1:2),INTENT(IN)::GAP
    INTEGER,INTENT(INOUT)::stateinfo
    REAL(RK_SOL),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(IN)::U
    REAL(RK_SOL),DIMENSION(S1:S2),INTENT(IN)::X
    REAL(RK_SOL),DIMENSION(1:INTER_Y),INTENT(IN)::Y
    REAL(RK_SOL),DIMENSION(:),ALLOCATABLE::UEXA,DIFF
    integer::i,j,k
    REAL(RK_SOL)::ErrN2_RED,ErrN2,Ninf_RED,Ninf
    CHARACTER(len=3)::rank,CAS
    !INITIALISATION DE UEXA SUIVANT LE CAS SOURCE
    ALLOCATE(UEXA(it1:itN))
          WRITE(rank,fmt='(1I3)')rang;WRITE(CAS,fmt='(1I3)')userchoice
          OPEN(10+rang,file='EXACTE_ERREUR/SOL_EXACTE_'//trim(adjustl(CAS))//'_ME'//trim(adjustl(rank))//'.dat')
    SELECT CASE(userchoice)
       CASE(1)!F1
          DO i=S1,S2
             DO j=GAP(i,1),GAP(i,2)
                k=j+(i-1)*INTER_Y
                UEXA(k)=X(i)*(1-X(i))*Y(j)*(1-Y(j))
                WRITE(10+rang,*)X(i),Y(j),UEXA(k)
             END DO
          END DO
       CASE(2)!F2
          DO i=S1,S2
             DO j=GAP(i,1),GAP(i,2)
                k=j+(i-1)*INTER_Y
                UEXA(k)=sin(X(i))+cos(Y(j))
                WRITE(10+rang,*)X(i),Y(j),UEXA(k)
             END DO
          END DO
       CASE DEFAULT!F3 PAS DE SOL EXACTE
          PRINT*,'PAS DE SOLUTION EXACTE POUR CE CAS (F3)'
    END SELECT
    CLOSE(10+rang)
   
       !ERREUR NORME_2:
       ALLOCATE(DIFF(it1:itN))
       DIFF=UEXA-U(it1:itN)
       ErrN2_RED=DOT_PRODUCT(DIFF,DIFF)
       CALL MPI_ALLREDUCE(ErrN2_RED,ErrN2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
       ErrN2=SQRT(ErrN2)
       !ERREUR NORME INFINIE
       Ninf_RED=MAXVAL(ABS(UEXA(it1:itN)-U(it1:itN)))
       CALL MPI_ALLREDUCE(Ninf_RED,Ninf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
       PRINT*,'ERREUR NORME 2:=',ErrN2,' ERREUR NORME INFINIE:=',Ninf
       OPEN(20+rang,file='EXACTE_ERREUR/ErrABSOLUE'//trim(adjustl(CAS))//'_ME'//trim(adjustl(rank))//'.dat')
       DO i=S1,S2
          DO j=GAP(i,1),GAP(i,2)
             k=j+(i-1)*INTER_Y
             WRITE(20+rang,*)X(i),Y(j),-DIFF(k)
          END DO
       END DO
       CLOSE(20+rang)
       DEALLOCATE(DIFF,UEXA)
   
    
  END SUBROUTINE UEXACT


  SUBROUTINE WR_U(U,it1,itN,INTER_X,INTER_Y,S1,S2,GAP,userchoice,rang,&
       IT_TEMPS,IT_MAX,Lx,Ly,X,Y)
    INTEGER,INTENT(IN)::userchoice,it1,itN,INTER_Y,INTER_X,S1,S2,rang,IT_TEMPS,IT_MAX
    INTEGER,DIMENSION(S1:S2,1:2),INTENT(IN)::GAP
    REAL(RK_SOL),INTENT(IN)::Lx,Ly
    REAL(RK_SOL),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(IN)::U
    REAL(RK_SOL),DIMENSION(S1:S2),INTENT(IN)::X
    REAL(RK_SOL),DIMENSION(1:INTER_Y),INTENT(IN)::Y
    integer::i,j,k
    CHARACTER(len=3)::rank,CAS
    CHARACTER(len=3)::Ntps
    WRITE(rank,fmt='(1I3)')rang;WRITE(CAS,fmt='(1I3)')userchoice;WRITE(Ntps,fmt='(1I3)')IT_TEMPS
    OPEN(10+rang,file='SOL_NUMERIQUE/U_ME'//trim(adjustl(rank))&
         //'_T'//trim(adjustl(Ntps)))
    SELECT CASE(userchoice)
       CASE(1)
          DO i=S1,S2
             DO j=GAP(i,1),GAP(i,2)
                k=j+(i-1)*INTER_Y
                IF(i-1==0)THEN
                   WRITE(10+rang,*)0.0d0,Y(j),0.0d0
                END IF
                IF(i+1==INTER_X+1)THEN
                   WRITE(10+rang,*)Lx,Y(j),0.0d0
                END IF
                IF(j-1==0)THEN
                   WRITE(10+rang,*)X(i),0.0d0,0.0d0
                END IF
                IF(j+1==INTER_Y+1)THEN
                   WRITE(10+rang,*)X(i),Ly,0.0d0
                END IF
                WRITE(10+rang,*)X(i),Y(j),U(k)
             END DO
          END DO
       CASE(2)
          DO i=S1,S2
             DO j=GAP(i,1),GAP(i,2)
                k=j+(i-1)*INTER_Y
                IF(i==1)THEN!BC OUEST (H)
                   WRITE(10+rang,*)0.0d0,Y(j),cos(Y(j))
                END IF
                IF(i==INTER_X)THEN!BC EST (H)
                   WRITE(10+rang,*)Lx,Y(j),cos(Y(j))+sin(Lx)
                END IF
                IF(j==1)THEN!BC SUD (G)
                   WRITE(10+rang,*)X(i),0.0d0,1.0d0+sin(X(i))
                END IF
                IF(j==INTER_Y)THEN!BC NORD (G)
                   WRITE(10+rang,*)X(i),Ly,cos(Ly)+sin(X(i))
                END IF
                WRITE(10+rang,*)X(i),Y(j),U(k)
             END DO
          END DO
       CASE(3)
          DO i=S1,S2
             DO j=GAP(i,1),GAP(i,2)
                k=j+(i-1)*INTER_Y
                IF(i-1==0)THEN!BC OUEST OU EST (H)
                   WRITE(10+rang,*)0.0d0,Y(j),1.0d0
                END IF
                IF(i+1==INTER_X+1)THEN!BC OUEST OU EST (H)
                   WRITE(10+rang,*)Lx,Y(j),1.0d0
                END IF
                IF(j-1==0.OR.j+1==INTER_Y+1)THEN!BC SUD OU NORD (G)
                   WRITE(10+rang,*)X(i),Y(j),0.0d0
                END IF
                IF(j+1==INTER_Y+1)THEN!BC SUD OU NORD (G)
                   WRITE(10+rang,*)X(i),Ly,0.0d0
                END IF
                WRITE(10+rang,*)X(i),Y(j),U(k)
             END DO
          END DO
    END SELECT
    CLOSE(10+rang)
  END SUBROUTINE WR_U


  SUBROUTINE VIZU_SOL_NUM(userchoice,DT,ITMAX,Nx,Ny)
    INTEGER,INTENT(INOUT)::userchoice,ITMAX,Nx,Ny
    REAL(RK_SOL),INTENT(IN)::DT
    CHARACTER(LEN=12)::VIZU
    INTEGER::i1,i2,i

    VIZU='VIZU_PLT.plt'

    PRINT*,'*****************************************************************************************'
    PRINT*,'************** VISUALISATION SOLUTION POUR LE CAS #',userchoice,' ***********************'
    PRINT*,'*****************************************************************************************'
    PRINT*,'************* SACHANT QUE DT=',DT,'ET ITMAX=',ITMAX,':'
    PRINT*,'*****************************************************************************************'
    PRINT*,'** ENTREZ LE NUMERO DE L''ITERATION A LAQUELLE VOUS VOULEZ COMMENCER LA VISUALISATION  **'
    READ*,i1
    PRINT*,'*****************************************************************************************'
    PRINT*,'*****************************************************************************************'
    PRINT*,'*** ENTREZ LE NUMERO DE L''ITERATION A LAQUELLE VOUS VOULEZ STOPPER LA VISUALISATION  ***'
    READ*,i2
    PRINT*,'*****************************************************************************************'
    PRINT*,'*****************************************************************************************'
    OPEN(10,file=VIZU)
    SELECT CASE(userchoice)
    CASE(1)
       WRITE(10,*)'set cbrange[0:1]'
    CASE(2)
       WRITE(10,*)'set cbrange[0.8:2]'
    CASE(3)
       WRITE(10,*)'set cbrange[0:1]' 
    END SELECT
    WRITE(10,*)'set dgrid3d,',Nx,',',Ny; WRITE(10,*)'set hidden3d'; WRITE(10,*)'set pm3d'
    WRITE(10,*)'do for [i=',i1,':',i2,']{splot ''SOL_NUMERIQUE/U.dat'' index i u 1:2:3 with pm3d at sb}'
    CLOSE(10)
    call system('gnuplot -p VIZU_PLT.plt')
  END SUBROUTINE VIZU_SOL_NUM


end module ModSOL
