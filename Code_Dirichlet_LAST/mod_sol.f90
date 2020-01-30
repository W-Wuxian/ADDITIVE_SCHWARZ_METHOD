module mod_sol
  use mpi
  use mod_parametres
  IMPLICIT NONE

  !*** MODULE ModSOL POUR DONNER LA SOLUTION EXACTE ET
  !*** POUR GERER L'ENREGISTREMENT DE LA SOLUTION DANS DES FICHIERS
  !*** SUBROUTINE POUR APPELLER UN SCRIPT GNUPLOT POUR LA VISUALISATION


CONTAINS

  SUBROUTINE UEXACT(U,userchoice)!userchoice==CT
    INTEGER,INTENT(INOUT)::userchoice
    REAL(PR),DIMENSION(it1:itN),INTENT(IN)::U
    REAL(PR),DIMENSION(:),ALLOCATABLE::UEXA,DIFF
    REAL(PR)::ErrN2_RED,ErrN2,Ninf_RED,Ninf
    CHARACTER(len=3)::CAS!rank,CAS
    !INITIALISATION DE UEXA SUIVANT LE CAS SOURCE
    ALLOCATE(UEXA(it1:itN))
    WRITE(CAS,fmt='(1I3)')userchoice
    OPEN(TAG,file='EXACTE_ERREUR/SOL_EXACTE_CAS'//trim(adjustl(CAS))//'_ME'//trim(adjustl(rank))//'.dat')
    SELECT CASE(userchoice)
    CASE(1)!F1
       DO i=S1,S2
          DO j=1,Ny_g 
             k=j+(i-1)*Ny_g
             UEXA(k)=X(i)*(1-X(i))*Y(j)*(1-Y(j))
             WRITE(TAG,*)X(i),Y(j),UEXA(k)
          END DO
       END DO
    CASE(2)!F2
       DO i=S1,S2
          DO j=1,Ny_g
             k=j+(i-1)*Ny_g
             UEXA(k)=sin(X(i))+cos(Y(j))
             WRITE(TAG,*)X(i),Y(j),UEXA(k)
          END DO
       END DO
    CASE DEFAULT!F3 PAS DE SOL EXACTE
       PRINT*,'PAS DE SOLUTION EXACTE POUR CE CAS (F3)'
    END SELECT
    CLOSE(TAG)

    !SUR LE DOMAINE [|S1;S2|]x[|1;Ny_g|]:
    !ERREUR NORME_2:
    ALLOCATE(DIFF(it1:itN))
    DIFF=UEXA-U(it1:itN)
    ErrN2_RED=DOT_PRODUCT(DIFF,DIFF)
    CALL MPI_ALLREDUCE(ErrN2_RED,ErrN2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
    ErrN2=SQRT(ErrN2)
    !ERREUR NORME INFINIE
    Ninf_RED=MAXVAL(ABS(UEXA(it1:itN)-U(it1:itN)))
    CALL MPI_ALLREDUCE(Ninf_RED,Ninf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
    PRINT*,'PAR DOMAINE: ','ERREUR NORME 2 :=',ErrN2,' ERREUR NORME INFINIE :=',Ninf
    OPEN(TAG,file='EXACTE_ERREUR/ErrABSOLUE_PAR_DOMAINE'//trim(adjustl(CAS))//'_ME'//trim(adjustl(rank))//'.dat')
    DO i=S1,S2
       DO j=1,Ny_g!GAP(i,1),GAP(i,2)
          k=j+(i-1)*Ny_g!Ny_g
          WRITE(TAG,*)X(i),Y(j),-DIFF(k)
       END DO
    END DO
    CLOSE(TAG)
    !SUR LE DOMAINE [|S1_old;S2_old|]x[|1;Ny_g|]:
    ErrN2_RED=0.0_PR; ErrN2=0.0_PR; Ninf_RED=0.0_PR; Ninf=0.0_PR
    !DIFF(it1_old:itN_old)=UEXA(it1_old:itN_old)-U(it1_old:itN_old)
    ErrN2_RED=DOT_PRODUCT(DIFF(it1_old:itN_old),DIFF(it1_old:itN_old))
    CALL MPI_ALLREDUCE(ErrN2_RED,ErrN2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
    ErrN2=SQRT(ErrN2)
    Ninf_RED=MAXVAL(ABS(UEXA(it1_old:itN_old)-U(it1_old:itN_old)))
    CALL MPI_ALLREDUCE(Ninf_RED,Ninf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
    PRINT*,'SUR [|S1old;S2old|]: ','ERREUR NORME 2 :=',ErrN2,' ERREUR NORME INFINIE :=',Ninf
    OPEN(TAG,file='EXACTE_ERREUR/ErrABSOLUE_oldDD'//trim(adjustl(CAS))//'_ME'//trim(adjustl(rank))//'.dat')
    DO i=S1_old,S2_old
       DO j=1,Ny_g!GAP(i,1),GAP(i,2)
          k=j+(i-1)*Ny_g!Ny_g
          WRITE(TAG,*)X(i),Y(j),-DIFF(k)
       END DO
    END DO
    CLOSE(TAG)
    DEALLOCATE(DIFF,UEXA)

  END SUBROUTINE UEXACT


  SUBROUTINE WR_U(U,IT_TEMPS)
    INTEGER,INTENT(IN)::IT_TEMPS
    REAL(PR),DIMENSION(it1_old:itN_old),INTENT(INOUT)::U
    CHARACTER(len=3)::CAS
    CHARACTER(len=3)::Ntps
    INTEGER::i,j
    WRITE(CAS,fmt='(1I3)')CT;WRITE(Ntps,fmt='(1I3)')IT_TEMPS
    OPEN(10+rang,file='SOL_NUMERIQUE/U_ME'//trim(adjustl(rank))&
         //'_T'//trim(adjustl(Ntps)))
    SELECT CASE(CT)
    CASE(1)
       DO i=S1_old,S2_old
          DO j=1,Ny_g
             k=j+(i-1)*Ny_g
             IF(i-1==0)THEN
                WRITE(10+rang,*)0.0_PR,Y(j),0.0_PR
             END IF
             IF(i+1==Nx_g+1)THEN
                WRITE(10+rang,*)Lx,Y(j),0.0_PR
             END IF
             IF(j-1==0)THEN
                WRITE(10+rang,*)X(i),0.0_PR,0.0_PR
             END IF
             IF(j+1==Ny_g+1)THEN
                WRITE(10+rang,*)X(i),Ly,0.0_PR
             END IF
             WRITE(10+rang,*)X(i),Y(j),U(k)
          END DO
       END DO
    CASE(2)
       DO i=S1_old,S2_old
          DO j=1,Ny_g
             k=j+(i-1)*Ny_g
             IF(i==1)THEN!BC OUEST (H)
                WRITE(10+rang,*)0.0_PR,Y(j),cos(Y(j))
             END IF
             IF(i==Nx_g)THEN!BC EST (H)
                WRITE(10+rang,*)Lx,Y(j),cos(Y(j))+sin(Lx)
             END IF
             IF(j==1)THEN!BC SUD (G)
                WRITE(10+rang,*)X(i),0.0_PR,1.0_PR+sin(X(i))
             END IF
             IF(j==Ny_g)THEN!BC NORD (G)
                WRITE(10+rang,*)X(i),Ly,cos(Ly)+sin(X(i))
             END IF
             WRITE(10+rang,*)X(i),Y(j),U(k)
          END DO
       END DO
    CASE(3)
       DO i=S1_old,S2_old
          DO j=1,Ny_g
             k=j+(i-1)*Ny_g
             IF(i-1==0)THEN!BC OUEST OU EST (H)
                WRITE(10+rang,*)0.0_PR,Y(j),1.0_PR
             END IF
             IF(i+1==Nx_g+1)THEN!BC OUEST OU EST (H)
                WRITE(10+rang,*)Lx,Y(j),1.0_PR
             END IF
             IF(j-1==0)THEN!BC SUD OU NORD (G)
                WRITE(10+rang,*)X(i),0.0_PR,0.0_PR
             END IF
             IF(j+1==Ny_g+1)THEN!BC SUD OU NORD (G)
                WRITE(10+rang,*)X(i),Ly,0.0_PR
             END IF
             WRITE(10+rang,*)X(i),Y(j),U(k)
          END DO
       END DO
    END SELECT
    CLOSE(10+rang)
  END SUBROUTINE WR_U


  SUBROUTINE VIZU_SOL_NUM(userchoice)!sequentiel call by proc 0 at the end
    INTEGER,INTENT(INOUT)::userchoice
    
    CHARACTER(LEN=12)::VIZU
    INTEGER::vi1,vi2

    VIZU='VIZU_PLT.plt'

    PRINT*,'*****************************************************************************************'
    PRINT*,'************** VISUALISATION SOLUTION POUR LE CAS #',userchoice,' ***********************'
    PRINT*,'*****************************************************************************************'
    PRINT*,'************* SACHANT QUE DT=',dt,'ET ITMAX=',Nt,':'
    PRINT*,'*****************************************************************************************'
    PRINT*,'** ENTREZ LE NUMERO DE L''ITERATION A LAQUELLE VOUS VOULEZ COMMENCER LA VISUALISATION  **'
    READ*,vi1
    PRINT*,'*****************************************************************************************'
    PRINT*,'*****************************************************************************************'
    PRINT*,'*** ENTREZ LE NUMERO DE L''ITERATION A LAQUELLE VOUS VOULEZ STOPPER LA VISUALISATION  ***'
    READ*,vi2
    PRINT*,'*****************************************************************************************'
    PRINT*,'*****************************************************************************************'
    OPEN(50,file=VIZU)
    SELECT CASE(userchoice)
    CASE(1)
       WRITE(50,*)'set cbrange[0:1]'
    CASE(2)
       WRITE(50,*)'set cbrange[0.8:2]'
    CASE(3)
       WRITE(50,*)'set cbrange[0:1]' 
    END SELECT
    WRITE(50,*)'set dgrid3d,',Nx_g+2,',',Ny_g+2; WRITE(50,*)'set hidden3d'; WRITE(50,*)'set pm3d'
    WRITE(50,*)'do for [i=',vi1,':',vi2,']{splot ''SOL_NUMERIQUE/U.dat'' index i u 1:2:3 with pm3d at sb}'
    CLOSE(50)
    call system('gnuplot -p VIZU_PLT.plt')
  END SUBROUTINE VIZU_SOL_NUM



  SUBROUTINE VIZU_SOL_NUM_ONLY_Nt(userchoice,MINU,MAXU)!sequentiel call by proc 0 at the end
    INTEGER,INTENT(INOUT) :: userchoice
    REAL(PR),INTENT(IN)   :: MINU, MAXU    
    CHARACTER(LEN=12)::VIZU
    

    VIZU='VIZU_PLT.plt'

    PRINT*,'*****************************************************************************************'
    PRINT*,'************** VISUALISATION SOLUTION POUR LE CAS #',userchoice,' ***********************'
    PRINT*,'*****************************************************************************************'
    PRINT*,'************* SACHANT QUE DT=',dt,'ET ITMAX=',Nt,':'
    PRINT*,'*****************************************************************************************'
    PRINT*,'** VISUALISATION AU TEMPS FINAL  **'
    PRINT*, (Nt)*dt,' secondes'
    PRINT*,'*****************************************************************************************'
    OPEN(50,file=VIZU)
    SELECT CASE(userchoice)
    CASE(1)
       WRITE(50,*)'set cbrange[',MINU,':',MAXU,']'
    CASE(2)
       WRITE(50,*)'set cbrange[',MINU,':',MAXU,']'
    CASE(3)
       WRITE(50,*)'set cbrange[',MINU,':',MAXU,']'
    END SELECT
    WRITE(50,*)'set dgrid3d,',Nx_g+2,',',Ny_g+2; WRITE(50,*)'set hidden3d'; WRITE(50,*)'set pm3d'
    WRITE(50,*)'splot ''SOL_NUMERIQUE/U.dat''  u 1:2:3 with pm3d at sb'
    CLOSE(50)
    call system('gnuplot -p VIZU_PLT.plt')
  END SUBROUTINE VIZU_SOL_NUM_ONLY_Nt


end module mod_sol
