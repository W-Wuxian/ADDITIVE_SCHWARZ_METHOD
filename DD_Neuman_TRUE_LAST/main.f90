!-----------------------------------
! MAIN : PROGRAMME PRINCIPAL
!-----------------------------------

Program main
!-------------modules-----------------!
  USE mpi
  use mod_parametres
  use mod_charge
  use mod_comm
  use mod_sol
  use mod_fonctions
  use mod_gradconj
  use mod_gmres
  use mod_matsys
!--------------------------------------!
  Implicit None
INTEGER :: DIMKRYLOV
INTEGER :: ILOOP
!VARIABLES POUR SAVE LE TEMPS DE CALCUL:
REAL(PR) :: TPS1, TPS2, MAX_TPS
REAL(PR) :: MINU, MAXU
REAL(PR) :: RED_MINU, RED_MAXU
REAL(PR) :: RED_SUMC, MAXC, TPSC1, TPSC2
REAL(PR) :: STA, STB, STM

!------------INIT REGION PARA-----------------!
!********* DEBUT REGION PARA:
CALL MPI_INIT(stateinfo) !INITIALISATION DU //LISME
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,stateinfo) !ON RECUPERE LES RANGS
!AU SEIN DE L'ENSEMBLE MPI_COMM_WORLD
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,Np,stateinfo)   !ET LE NOMBRE DU PROCESSUS
!---------------------------------------------!

!******** READ USERCHOICE FOR CT AND SHARE IT WITH OTHER PROCESS
IF(rang==0)THEN
   CALL SYSTEM('make cleanREP')
   Print*,"______RESOLUTION DE L'ÉQUATION : DD SCHWARZ ADDITIVES //______"
   !---initialisation des variables cas test-------------
   Print*,"Donnez le cas test (1,2 ou 3):"
   Read*,CT
   PRINT*," POUR VISUALISATION ENTREZ 1, SINON 0"
   read*,sysmove
END IF
CALL MPI_BCAST(CT,1,MPI_INTEGER,0,MPI_COMM_WORLD,stateinfo)
CALL MPI_BCAST(sysmove,1,MPI_INTEGER,0,MPI_COMM_WORLD,stateinfo)

TPS1=MPI_WTIME()

!******** COMPUTE SOME MATRIX_G DATA & GEO TIME DATA
nnz_g=5*Na_g-2*Nx_g-2*Ny_g    !nb d'elt non nul dans la matrice A_g
IF(rang==0)THEN
   print*,"Nombre éléments non nuls du domaine A_g non DD ADDITIVE : nnz_g=",nnz_g
   print*,'Nombre de ligne Nx_g, colonne Ny_g de A_g',Nx_g,Ny_g
   Print*,"Lx=",Lx
   Print*,"Ly=",Ly
   Print*,"D=",D
   Print*,"   "
   Print*,"dx=",dx
   Print*,"dy=",dy
   Print*,"dt",dt
   Print*,"Tfinal=",Tf,"secondes"
   Print*,"   "
   print*,"alpha=",alpha
   print*,"beta=",beta
   print*,"gamma=",gamma
   Print*,"   "
END IF

!******* SET UP ONE TAG FILE NAMED TAG TO i/o:
TAG=10+rang
WRITE(rank,fmt='(1I3)')rang !WATCH OUT IF rang>999

!******* COMPUTE SIZE OF LOCAL DD OVER LOCAL MATRIX A_l SIZE:
OPEN(TAG,file='CHARGE/CHARGE'//trim(adjustl(rank))//'.dat')
CALL MODE1(rang, Np, S1, S2, it1, itN)!-->DISTRIBUTION DE LA CHARGE PAR PROCESSUS
S1_old = S1; S2_old = S2; it1_old = it1; itN_old = itN
WRITE(TAG,*)'AVANT OVERLAP', S1, S2, it1, itN
CALL CHARGE_OVERLAP(Np, rang, S1, S2, it1, itN)!-->MODIFICATION DE LA CHARGE POUR DIRICHLET
WRITE(TAG,*)'APRES OVERLAP DIRICHLET', S1, S2, it1, itN
CALL CHARGE_NEUMAN(it1, itN, rang, Np, S1, S2)!-->MODIFICATION DE LA CHARGE CALCULEE PAR (CHARGE_OVERLAP) POUR NEUMANN
WRITE(TAG,*)'OVERLAP FINAL NEUMAN', S1, S2, it1, itN

Nx_l = S2 - S1 + 1; Ny_l = Ny_g; Na_l = Nx_l * Ny_l
nnz_l = 5 * Na_l - 2 * Nx_l - 2 * Ny_l
crtl_nnz_l = (Nx_l*Ny_l) + 2 * (Nx_l) * (Ny_l - 1) + 2 * (Nx_l - 1) * (Ny_l)
WRITE(TAG,*)'Nxl,Ny_l,Na_l', Nx_l, Ny_l, Na_l
WRITE(TAG,*)'nnz_l,crtl_nnz_l', nnz_l,crtl_nnz_l
CLOSE(TAG)
!CALCUL DES DIMMENSIONS DES MATRICES:
D_rows = Ny_l; D_nnz = D_rows+2*(D_rows-1)
C_rows = Ny_l; C_nnz = C_rows
crtl_L1_nnz = D_nnz+C_nnz
crtl_L2_nnz = (Nx_l-2)*(D_nnz+2*C_nnz)
crtl_L3_nnz = D_nnz+C_nnz
sum_crtl_L_nnz = crtl_L1_nnz+crtl_L2_nnz+crtl_L3_nnz
crtl_L1_nnz_Nxl_EQV_1 = D_nnz

!test pour la verification du nombre de non zero des matrices locales:
IF(Nx_l > 1)THEN
   IF(crtl_nnz_l/=nnz_l .AND. nnz_l/=sum_crtl_L_nnz)THEN
      PRINT*,'ERROR,RANG=',rang,'Local Matrix A_l.rows,A_l.cols',Nx_l,Ny_l,nnz_l,&
           &'A_l_nnz,crtl_nnz_l,sum_crtl_L_nnz',nnz_l,crtl_nnz_l,sum_crtl_L_nnz
      GOTO 9999
   END IF
ELSE IF(Nx_l == 1)THEN
   IF(crtl_nnz_l/=nnz_l .AND. nnz_l/=crtl_L1_nnz_Nxl_EQV_1)THEN
      PRINT*,'ERROR,RANG=',rang,'Local Matrix A_l.rows,A_l.cols',Nx_l,Ny_l,nnz_l,&
           &'A_l_nnz,crtl_nnz_l,crtl_L1_nnz_Nxl_EQV_1',nnz_l,crtl_nnz_l,crtl_L1_nnz_Nxl_EQV_1
      GOTO 9999
   END IF
END IF
!******* COMPUTE CARTESIAN GRID OVER PROCESS i.e CHARGE APPLY TO X,Y,T
IF(Np/=1)THEN
   CALL param_para(rang,Np,S1,S2,Ny_l,X,Y,T)
ELSE
   CALL param(X,Y,T)
   PRINT*,'EN SEQUENTIEL'
END IF

!****** ALLOCATE Uo ET U
ALLOCATE(U(it1:itN),Uo(it1:itN))
Uo = 1._PR      !CI pour le gradient conjugué
U  = 1._PR       !vecteur solution initialisé
!SAVE SOL :
call WR_U(U,0)
CALL matsys_v2(nnz_l,Nx_l,Ny_l,AA,IA,JA)!-->on construit A au format COO (AA,IA,JA)

OPEN(TAG,file='MAT_VIZU/matrice_A_'//trim(adjustl(rank))//'.dat')!-->POUR VISUALISER LES MATRICES
DO i=1,nnz_l
  WRITE(TAG,*)AA(i),IA(i),JA(i)
END DO
CLOSE(TAG)

ALLOCATE(F(it1:itN))!-->SECOND MEMBRE
ALLOCATE(UG(LBUG:UBUG),UD(LBUD:UBUD))
UG = 0.0_PR; UD = 0.0_PR !-->INITITALISATION
!---boucle principale (résolution systeme) et ecriture a chaque itération---
Print*,"----------DÉBUT DU CALCUL DE LA SOLUTION------------"

Pivot=(Np-mod(Np,2))/2 !-->!PIVOT DE COMMUNICATION:
RED_SUMC = 0.0_PR

DIMKRYLOV = 5!Na_l-5
IF(Np>1)THEN
   boucle_temps:Do ILOOP=1,Nt   !-->BOUCLE EN TEMPS
      err = 1.0_PR; step = 1; Norm2 = 0.0_PR
      !STA = MPI_WTIME()
      boucle_CV_SWHARZ:Do while(err>err_LIM .OR. step<3)!-->BOUCLE POUR CONVERGENCE DE SCHWARZ
         step = step +1
         Norm1 = Norm2
         !Norm2 = 0.0_PR
         call vectsource_FULL(CT,U,UG,UD,S1,S2,X,Y,T(ILOOP),F)!-->SECOND MEMBRE F=Uo+dt*SOURCE_TERME+CL
         !call MOD_gradconj_SPARSE(AA,IA,JA,f,Uo,U,res,k,it1,itN)
         !call MOD_jacobi_SPARSE(AA,IA,JA,F,Uo,U,res,k)
         call BICGSTAB(it1, itN, AA, IA, JA, U, F, 0.0000000001_PR)! solver le plus adapté au PB
         !call GMres_COO(AA,IA,JA,F,Uo,U,res,k,Na_l,it1,itN,DIMKRYLOV)! mal interfacé

         TPSC1 = MPI_WTIME()
         !CALL COMM_PTOP_2WAY(U,UG,UD)!-->UNE PETITE REDUCTION DU TEMPS DE COMMUNICATION(NP>=5)
         call COMM_PTOP(U,UG,UD)!-->COMMUNICATION DES FRONTIERES "IMMERGEE"
         TPSC2 = MPI_WTIME()
         RED_SUMC = RED_SUMC + (TPSC2-TPSC1)

         Uo=U !-->SWAP
         !POUR LE CALCUL DE L'ERREUR DE SCHWARZ(C1,Nomr1,Norm2,err):
         if(rang==0)then
            C1 = DOT_PRODUCT(U(itN-Ny_l+1:itN),U(itN-Ny_l+1:itN))
         else if(rang==Np-1)then
            C1 = DOT_PRODUCT(U(it1:it1+Ny_l-1),U(it1:it1+Ny_l-1))
         else
            C1 = DOT_PRODUCT(U(it1:it1+Ny_l-1),U(it1:it1+Ny_l-1)) +DOT_PRODUCT(U(itN-Ny_l+1:itN),U(itN-Ny_l+1:itN))
         end if
         !C1 = DOT_PRODUCT(U(it1:itN),U(it1:itN))!-->AUTRE METHODE POUR LE CALCUL DE C1
         Call MPI_ALLreduce(C1,Norm2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
         Norm2 = SQRT(Norm2)
         err = abs(Norm1-Norm2)
         !PRINT*,'err=',err,'step',step,'Norm2',Norm2,'Norm1',Norm1
      END DO boucle_CV_SWHARZ
      !STB = MPI_WTIME()
      if(mod(ILOOP,10)==0)then
         Print*,'IT TEMPS',ILOOP
      end if
   End DO boucle_temps
   !CALL MPI_ALLREDUCE(STB-STA,STM,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,stateinfo)
   !PRINT*,'STM',STM,'STEP',STEP
ELSE IF(Np==1)THEN
   bboucle_temps:Do ILOOP=1,Nt   
         call vectsource_FULL(CT,U,UG,UD,S1,S2,X,Y,T(ILOOP),F)
         !call MOD_gradconj_SPARSE(AA,IA,JA,f,Uo,U,res,k,it1,itN)
         !call MOD_jacobi_SPARSE(AA,IA,JA,F,Uo,U,res,k)
         call BICGSTAB(it1, itN, AA, IA, JA, U, F, 0.000000000001_PR)
         !call GMres_COO(AA,IA,JA,F,Uo,U,res,k,Na_l,it1,itN,DIMKRYLOV)

         Uo=U !SWAP
         if(mod(ILOOP,10)==0)then
            Print*,'IT TEMPS',ILOOP
         end if
   End DO bboucle_temps
END IF

!----------------------------------------------
!SAVE SOL :
call WR_U(U,Nt)
DEALLOCATE(UG,UD)

!WRITE SOL EXACTE:
IF(CT<3)THEN
   call  UEXACT(U,CT)
END IF

RED_MINU = minval(U)
RED_MAXU = maxval(U)
CALL MPI_ALLREDUCE(RED_MINU,MINU,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,stateinfo)
CALL MPI_ALLREDUCE(RED_MAXU,MAXU,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
if(Np>1)then
   CALL MPI_ALLREDUCE(RED_SUMC,MAXC,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
end if


DEALLOCATE(U,Uo,X,Y,T,AA,IA,JA,F)

TPS2=MPI_WTIME()
CALL MPI_ALLREDUCE(TPS2-TPS1,MAX_TPS,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
IF(rang==0)THEN
   PRINT*,'TEMPS CALCUL=',MAX_TPS,' TEMPS COMM=',MAXC
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
IF(crtl_nnz_l == nnz_l)THEN
   GOTO 6666
END IF
9999 PRINT*,rang,'GOTO9999crtl_nnz_l /= nnz_l ERROR',crtl_nnz_l,nnz_l
6666 PRINT*,rang,'GOTO6666NEXT INSTRUCTION IS MPI_FINALIZE',crtl_nnz_l,nnz_l

!---ecriture solution finale--------------------------
!call CONCAT_VIZU(sysmove)
call CONCAT_VIZU_ONLY_Nt(sysmove,MINU,MAXU)

CALL MPI_FINALIZE(stateinfo)





CONTAINS

  Subroutine CONCAT_VIZU(sysmove)
    Integer, Intent(INOUT):: sysmove
    CHARACTER(LEN=3):: NNT
    CHARACTER(LEN=1)::ESPACE
    CHARACTER(LEN=20)::NFILE,OUT_FILE
    INTEGER::i
    !VIZUALISATION SEQUENTIELLE: ON CONCATENE POUR UNE ITERATION EN TEMPS DONNEE, LES FICHIERS DE RANG DIFFERENTS DANS UN FICHIER
    IF(rang==0.AND.sysmove==1)THEN
       ESPACE=' '
       DO i=Nt,Nt!0,Nt
          WRITE(NNT,fmt='(1I3)')i;
          NFILE='U_ME*'//'_T'//trim(adjustl(NNT));OUT_FILE='UT_'//trim(adjustl(NNT))//'.dat'
          !ON CONCATENE POUR UNE ITERATION EN TEMPS DONNEE, LES FICHIERS DE RANG DIFFERENTS DANS UN FICHIER
          CALL system('cat SOL_NUMERIQUE/'//NFILE//' >> SOL_NUMERIQUE/'//OUT_FILE)
          CALL system('rm SOL_NUMERIQUE/'//NFILE)
          OPEN(20,file='SOL_NUMERIQUE/'//OUT_FILE,position='append')
          SELECT CASE(CT)!ICI IL FAUT RAJOUTER LES POINTS AUX COINS DU DOMAINES
          CASE(1)
             WRITE(20,*)0.0d0,0.0d0,0.0d0;WRITE(20,*)Lx,0.0d0,0.0d0
             WRITE(20,*)0.0d0,Ly,0.0d0;WRITE(20,*)Lx,Ly,0.0d0
          CASE(2)
             WRITE(20,*)0.0d0,0.0d0,1.0d0;WRITE(20,*)Lx,0.0d0,1.0d0+sin(Lx)
             WRITE(20,*)0.0d0,Ly,cos(Ly);WRITE(20,*)Lx,Ly,cos(Ly)+sin(Lx)
          CASE(3)
             WRITE(20,*)0.0d0,0.0d0,0.5d0;WRITE(20,*)Lx,0.0d0,0.5d0!TEMPERATURES IMPOSEES DIFFERENTES, ON MOYENNE:((G+H)/2)
             WRITE(20,*)0.0d0,Ly,0.5d0;WRITE(20,*)Lx,Ly,0.5d0
          END SELECT
          WRITE(20,fmt='(1A1)')ESPACE
          WRITE(20,fmt='(1A1)')ESPACE
          CLOSE(20)
          !ON CONCATENE POUR UNE ITERATION EN TEMPS DONNEE, LES FICHIERS DE TEMPS DIFFERENTS DANS UN FICHIER
          CALL system('cat SOL_NUMERIQUE/'//OUT_FILE//' >> SOL_NUMERIQUE/U.dat')
          CALL system('rm SOL_NUMERIQUE/'//OUT_FILE)
       END DO
       !A LA FIN, IL RESTE UN FICHIER OU LA SOLUTION EN TEMPS EST SOUS FORME D'INDEX GNUPLOT
       !ON APPELLE LA SUBROUTIEN DE VISUALISATION:
       CALL VIZU_SOL_NUM(CT)
    END IF
  End Subroutine CONCAT_VIZU



  Subroutine CONCAT_VIZU_ONLY_Nt(sysmove,MINU,MAXU)
    Integer, Intent(INOUT):: sysmove
    REAL(PR), Intent(IN)              :: MINU, MAXU
    CHARACTER(LEN=3):: NNT
    CHARACTER(LEN=1)::ESPACE
    CHARACTER(LEN=20)::NFILE,OUT_FILE
    INTEGER::i
    !VIZUALISATION SEQUENTIELLE: ON CONCATENE POUR UNE ITERATION EN TEMPS DONNEE, LES FICHIERS DE RANG DIFFERENTS DANS UN FICHIER
    IF(rang==0.AND.sysmove==1)THEN
       ESPACE=' '
       DO i=Nt,Nt
          WRITE(NNT,fmt='(1I3)')i;     
          NFILE='U_ME*'//'_T'//trim(adjustl(NNT));OUT_FILE='UT_'//trim(adjustl(NNT))//'.dat'
          !ON CONCATENE POUR UNE ITERATION EN TEMPS DONNEE, LES FICHIERS DE RANG DIFFERENTS DANS UN FICHIER
          CALL system('cat SOL_NUMERIQUE/'//NFILE//' >> SOL_NUMERIQUE/'//OUT_FILE)
          CALL system('rm SOL_NUMERIQUE/'//NFILE)
          OPEN(20,file='SOL_NUMERIQUE/'//OUT_FILE,position='append')
          SELECT CASE(CT)!ICI IL FAUT RAJOUTER LES POINTS AUX COINS DU DOMAINES
          CASE(1)
             WRITE(20,*)0.0d0,0.0d0,0.0d0;WRITE(20,*)Lx,0.0d0,0.0d0
             WRITE(20,*)0.0d0,Ly,0.0d0;WRITE(20,*)Lx,Ly,0.0d0
          CASE(2)
             WRITE(20,*)0.0d0,0.0d0,1.0d0;WRITE(20,*)Lx,0.0d0,1.0d0+sin(Lx)
             WRITE(20,*)0.0d0,Ly,cos(Ly);WRITE(20,*)Lx,Ly,cos(Ly)+sin(Lx)
          CASE(3)
             WRITE(20,*)0.0d0,0.0d0,0.5d0;WRITE(20,*)Lx,0.0d0,0.5d0!TEMPERATURES IMPOSEES DIFFERENTES, ON MOYENNE:((G+H)/2)
             WRITE(20,*)0.0d0,Ly,0.5d0;WRITE(20,*)Lx,Ly,0.5d0
          END SELECT
          WRITE(20,fmt='(1A1)')ESPACE
          WRITE(20,fmt='(1A1)')ESPACE
          CLOSE(20)
          !ON CONCATENE POUR UNE ITERATION EN TEMPS DONNEE, LES FICHIERS DE TEMPS DIFFERENTS DANS UN FICHIER
          CALL system('cat SOL_NUMERIQUE/'//OUT_FILE//' >> SOL_NUMERIQUE/U.dat')
          CALL system('rm SOL_NUMERIQUE/'//OUT_FILE)
       END DO
       !A LA FIN, IL RESTE UN FICHIER OU LA SOLUTION EN TEMPS EST SOUS FORME D'INDEX GNUPLOT
       !ON APPELLE LA SUBROUTIEN DE VISUALISATION:
       !CALL VIZU_SOL_NUM(CT)
       CALL VIZU_SOL_NUM_ONLY_Nt(CT,MINU,MAXU)
    END IF
  End Subroutine CONCAT_VIZU_ONLY_Nt

End Program main
