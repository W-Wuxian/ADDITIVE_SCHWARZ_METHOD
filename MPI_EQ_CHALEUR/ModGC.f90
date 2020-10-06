module ModGC
  use mpi

  implicit none
  INTEGER,PARAMETER::RK_GC=8
  !ModGC CONTIENT: 
  !SUBROUTINE: GRADIENT CONJUGUE, 2 TYPES DE COMMUNICATION VECTEUR, PRODUIT MATRICE/VECTEUR

contains

  !GRADIENT CONJUGUE:
  SUBROUTINE GC(A,AP,AL,U,F,F_VAR,userchoice,DT,TOL,KMAX,it1,itN&
       ,INTER_Y,INTER_X,S1,S2,GAP,Np,rang,stateinfo,status,mode,r,z,p,&
       Tpmv,Tcomm,Tall,PIVOT,NUM_tps)
    REAL(RK_GC),INTENT(IN)::A,AP,AL
    INTEGER,INTENT(IN)::it1,itN,INTER_Y,INTER_X,S1,S2,Np,rang,userchoice,KMAX,mode,PIVOT,NUM_tps
    INTEGER,INTENT(IN),DIMENSION(S1:S2,1:2)::GAP
    REAL(RK_GC),INTENT(IN)::DT,TOL,F_VAR
    REAL(RK_GC),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(INOUT)::U
    REAL(RK_GC),DIMENSION(it1:itN),INTENT(IN)::F
    INTEGER,DIMENSION(MPI_STATUS_SIZE),INTENT(inout)::status
    INTEGER,INTENT(INOUT)::stateinfo
    REAL(RK_GC)::beta,beta_RED,alpha,beta_old,DOTzp_RED,DOTzp,gamma
    REAL(RK_GC),DIMENSION(it1:itN),INTENT(INOUT)::r,z
    REAL(RK_GC),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(INOUT)::p
    REAL(RK_GC),INTENT(INOUT)::Tpmv,Tcomm,Tall!TEMPS (SUM) produit matrice/vecteur,communication vecteur, MPI_ALLREDUCE
    REAL(RK_GC)::Tpmv_RED,Tcomm_RED,Tall_RED!VARIABLES DE REDUCTION
    REAL(RK_GC)::Tp1,Tp2,Tc1,Tc2,Ta1,Ta2!TEMPS
    !REAL(RK_GC),DIMENSION(:),ALLOCATABLE::r,p,z
    INTEGER::it_GC
    !ALLOCATE(r(it1:itN),z(it1:itN),p(it1-INTER_Y:itN+INTER_Y))
    !INITIALISATION:
    it_GC=0
    r=0.0d0;z=0.0d0
    !INITIALISATION TEMPS:
    Tp1=0.0d0;Tp2=0.0d0;Ta1=0.0d0;Ta2=0.0d0;Tc1=0.0d0;Tc2=0.0d0
    Tpmv_RED=0.0d0;Tcomm_RED=0.0d0;Tall_RED=0.0d0
    Tpmv=0.0d0;Tcomm=0.0d0;Tall=0.0d0
    !UN PRODUIT MATRICE/VECTEUR NECESSITE UNE COMM DE PARTIE DU VECTEUR ENTRE PROCS:
    Tc1=MPI_WTIME()
    IF(Np>1)THEN!DECOMMENTER LE MODE DE COMMUNICATION QUE VOUS VOULEZ, COMMENTER L'AUTRE
       CALL COMM_VECT(status,stateinfo,rang,Np,U,it1,itN,INTER_Y,INTER_X,mode)
       !CALL COMM_VECT_PIVOT(status,stateinfo,rang,Np,U,it1,itN,INTER_Y,INTER_X,mode,PIVOT)
    END IF
    Tc2=MPI_WTIME()
    Tp1=MPI_WTIME()
    CALL PMV(A,AP,AL,INTER_Y,INTER_X,it1,itN,r,U,S1,S2,GAP)
    Tp2=MPI_WTIME()
    !ON NE SOUHAITE PAS RECALCULER F3, NI STOCKER b=U+DT*F
    SELECT CASE(userchoice)
    CASE(3)!POUR F=F3, DANS CE CAS ON RAJOUTE LA VARIABLE EN TEMPS POUR LA FORCE F_VAR=cos(Pi*t/2)
       r=-r+U(it1:itN)+DT*F_VAR*F!.EQV. r=b-A.U avec b=U(it1:itN)+DT*F_VAR*F
    CASE DEFAULT!POUR F=F1 OU F=F2
       r=-r+U(it1:itN)+DT*F!.EQV. r=b-A.U avec b=U(it1:itN)+DT*F
    END SELECT
    
    beta_RED=DOT_PRODUCT(r,r)!NORME2 AU CARRE DU RESIDU (PAR MORCEAU)
    Ta1=MPI_WTIME()
    CALL MPI_ALLREDUCE(beta_RED,beta,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)!reduction
    Ta2=MPI_WTIME()
    
    p(it1-INTER_Y:it1-1)=0.0d0;p(it1:itN)=r(it1:itN);p(itN+1:itN+INTER_Y)=0.0d0!p la direction de descente du GC
    
    !TEMPS
    Tpmv_RED=Tp2-Tp1;Tcomm_RED=Tc2-Tc1;Tall_RED=Ta2-Ta1;
    DO WHILE(SQRT(beta)>TOL.AND.(it_GC<KMAX))!BOUCLE DE CONVERGENCE
       !z=Ap, OBLIGATION DE FAIRE UNE COMM SUR p:
       Tc1=MPI_WTIME()
       IF(Np>1)THEN!DECOMMENTER LE MODE DE COMMUNICATION QUE VOUS VOULEZ, COMMENTER L'AUTRE
          CALL COMM_VECT(status,stateinfo,rang,Np,p,it1,itN,INTER_Y,INTER_X,mode)
          !CALL COMM_VECT_PIVOT(status,stateinfo,rang,Np,p,it1,itN,INTER_Y,INTER_X,mode,PIVOT)
       END IF
       Tc2=MPI_WTIME()
       Tcomm_RED=Tcomm_RED+Tc2-Tc1
       Tp1=MPI_WTIME()
       CALL PMV(A,AP,AL,INTER_Y,INTER_X,it1,itN,z,p,S1,S2,GAP)!PRODUIT MATRICE VECTEUR
       Tp2=MPI_WTIME()
       Tpmv_RED=Tpmv_RED+Tp2-Tp1
       !REDUCTION SUR DOT_PRODUCT(z,p):
       DOTzp_RED=DOT_PRODUCT(z(it1:itN),p(it1:itN))
       Ta1=MPI_WTIME()
       CALL MPI_ALLREDUCE(DOTzp_RED,DOTzp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
       Ta2=MPI_WTIME()
       Tall_RED=Tall_RED+Ta2-Ta1
       !CALCUL DE alpha le pas de descente:
       alpha=beta/DOTzp

       !ACTUALISATION DE U ET de r:
       U(it1:itN)=U(it1:itN)+alpha*p(it1:itN)
       r=r-alpha*z

       !SAVE beta dans beta_old:
       beta_old=beta
       
       !ACTUALISATION DE BETA ET REDUCTION:
       beta_RED=DOT_PRODUCT(r,r)
       Ta1=MPI_WTIME()
       CALL MPI_ALLREDUCE(beta_RED,beta,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
       Ta2=MPI_WTIME()
       Tall_RED=Tall_RED+Ta2-Ta1
       !CALCUL de gamma:
       gamma=beta/beta_old
       
       !ACTUALISATION DE LA DIRECTION DE DESCENTE p:
       p(it1:itN)=r(it1:itN)+gamma*p(it1:itN)

       !ACTUALISATION DE k:
       it_GC=it_GC+1
       
    END DO
    !SORTIE DES MAX(SUM TEMPS):
    CALL MPI_ALLREDUCE(Tcomm_RED,Tcomm,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
    CALL MPI_ALLREDUCE(Tall_RED,Tall,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
    CALL MPI_ALLREDUCE(Tpmv_RED,Tpmv,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,stateinfo)
    
    PRINT*,'it temps',NUM_tps,'it_GC',it_GC,'SQRT(beta)',SQRT(beta),'TOL',TOL,'RANG',rang
    
   
  END SUBROUTINE GC



  !POUR COMMUNIQUER DES PARTIES D'UN VECTEUR:
  SUBROUTINE COMM_VECT(status,stateinfo,rang,NBR_P,ES,it1,itN,INTER_Y,INTER_X,mode)
    INTEGER,DIMENSION(MPI_STATUS_SIZE),INTENT(inout)::status
    INTEGER,INTENT(INOUT)::stateinfo
    INTEGER,INTENT(IN)::rang,NBR_P,it1,itN,INTER_Y,INTER_X,mode
    REAL(RK_GC),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(INOUT)::ES!LE VECTEUR A COMMUNIQUER
    INTEGER::NP,NPP,T1M,T1MM,REQUEST,A,B
    NP=itN+1;NPP=itN+INTER_Y;T1M=it1-1;T1MM=it1-INTER_Y
    A=itN-INTER_Y+1;B=it1+INTER_Y-1
       IF(NBR_P<=INTER_X)THEN
          IF(rang==0)THEN
             CALL MPI_SEND(ES(A:itN),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
             CALL MPI_RECV(ES(NP:NPP),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo) 
             !SEND PARTIE NORD DE ES AU RANG+1 ET RECV PARTIE SUD DE ES DU RANG+1

          ELSE IF(rang>0.AND.rang<NBR_P-1)THEN
             CALL MPI_RECV(ES(T1MM:T1M),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
             CALL MPI_SEND(ES(it1:B),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
             !RECV PARTIE NORD DE ES DU RANG-1 ET SEND PARTIE SUD DE ES AU RANG-1 
             CALL MPI_SEND(ES(A:itN),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
             CALL MPI_RECV(ES(NP:NPP),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
             !SEND PARTIE NORD DE ES AU RANG+1 ET RECV PARTIE SUD DE ES DU RANG+1                

          ELSE IF(rang==NBR_P-1)THEN
             CALL MPI_RECV(ES(T1MM:T1M),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
             CALL MPI_SEND(ES(it1:B),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
             !RECV PARTIE NORD DE ES DU RANG-1 ET SEND PARTIE SUD DE ES AU RANG-1
          END IF
          
       ELSE
          CALL COMM_NP_SUP(NBR_P,rang,it1,itN,INTER_Y,ES,stateinfo,status)!SI Np>INTER_X
       END IF
    
  END SUBROUTINE COMM_VECT

  ! POUR LE MODE 2 DANS LE CAS Np>Nx-2, un processus peut avoir a send et recv
  ! pour plus de 1 processus de rang sup ou de rang inf:
  ! SANS TRI obligation de faire comm point à point avec remonté et descente pour transmettre l'information 
  subroutine COMM_NP_SUP(Np,rang,it1,itN,INTER,ENTER,stateinfo,status)
    integer::k,wr!i,j,k,liste
    integer,intent(in)::Np,rang,it1,itN,INTER!INTER=INTER_Y

    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    real(RK_GC),dimension(it1-INTER:itN+INTER),intent(inout)::ENTER

    INTEGER::A,B,C,D,E,F
    A=itN-INTER+1;B=it1-INTER;C=it1-1;D=it1+INTER-1;E=itN+1;F=itN+INTER
    k=0
    wr=0
    
    !** BOUCLE COMMUNICATION POINT A POINT: RANG SEND PARTIE NORD DE ES A RANG+1 ET RANG+1 RECV PARTIE NORD DE ES  DE RANG
    do wr=0,Np-1
       if(wr==rang)then
          if(rang<Np-1)then
             call MPI_SEND(ENTER(A:itN),INTER,MPI_DOUBLE_PRECISION,rang+1,21,MPI_COMM_WORLD,stateinfo)
             !call MPI_SEND(ENTER(itN-INTER+1:itN),INTER,MPI_DOUBLE_PRECISION,rang+1,21,MPI_COMM_WORLD,stateinfo)
          end if
       else if(rang==wr+1)then
          call MPI_RECV(ENTER(B:C),INTER,MPI_DOUBLE_PRECISION,rang-1,21,MPI_COMM_WORLD,status,stateinfo)
          !call MPI_RECV(ENTER(it1-INTER:it1-1),INTER,MPI_DOUBLE_PRECISION,rang-1,21,MPI_COMM_WORLD,status,stateinfo)
       end if
    end do
    !** BOUCLE COMMUNICATION POINT A POINT: SENS INVERSE
    do wr=Np-1,0,-1
       if(wr==rang)then
          if(rang>0)then
             call MPI_SEND(ENTER(it1:D),INTER,MPI_DOUBLE_PRECISION,rang-1,31,MPI_COMM_WORLD,stateinfo)
             !call MPI_SEND(ENTER(it1:it1+INTER-1),INTER,MPI_DOUBLE_PRECISION,rang-1,31,MPI_COMM_WORLD,stateinfo)
          end if
       else if(rang==wr-1)then
          call MPI_RECV(ENTER(E:F),INTER,MPI_DOUBLE_PRECISION,rang+1,31,MPI_COMM_WORLD,status,stateinfo)
          !call MPI_RECV(ENTER(itN+1:itN+INTER),INTER,MPI_DOUBLE_PRECISION,rang+1,31,MPI_COMM_WORLD,status,stateinfo)
       end if
    end do

  end subroutine COMM_NP_SUP


!PRODUIT MATRICE VECTEUR PAR MORCEAU, ON NE STOCKE QUE LES TROIS COEFFICIANTS NON NULS DE LA MATRICE
!A,AP,AL, ETANT DONNE QUE LA MATRICE EST PENTADIAGONALE SDP, A COEFFICIANTS CONSTANTS PAR DIAGONALE
  SUBROUTINE PMV(A,AP,AL,INTER_Y,INTER_X,it1,itN,E,V,S1,S2,GAP)
    REAL(RK_GC),INTENT(IN)::A,AP,AL!COEF MATRICE
    INTEGER,INTENT(IN)::INTER_Y,INTER_X,it1,itN,S1,S2
    INTEGER,INTENT(IN),DIMENSION(S1:S2,1:2)::GAP
    real(RK_GC),dimension(it1:itN),intent(out)::E !E=A.V
    real(RK_GC),dimension(it1-INTER_Y:itN+INTER_Y),intent(out)::V
    INTEGER::i,j,k
    !BC=CONDITIONS LIMITES
    DO i=S1,S2
       IF(i==1)THEN!BC
          DO j=GAP(1,1),GAP(1,2)
             IF(j==1)THEN!BC
                E(1)=A*V(1)-AP*V(2)-AL*V(1+INTER_Y)
             ELSE IF(j<INTER_Y)THEN!BC
                E(j)=-AP*V(j-1)+A*V(j)-AP*V(j+1)-AL*V(j+INTER_Y)
             ELSE!BC
                E(j)=-AP*V(INTER_Y-1)+A*V(INTER_Y)-AL*V(2*INTER_Y)
             END IF
          END DO
       ELSE IF(i>1.AND.i<INTER_X)THEN
          DO j=GAP(i,1),GAP(i,2)
             k=j+(i-1)*INTER_Y
             IF(j==1)THEN!BC
                E(k)=-AL*V(k-INTER_Y)+A*V(k)-AP*V(k+1)-AL*V(k+INTER_Y)
             ELSE IF(j<INTER_Y)THEN!PAS DE BC 
                E(k)=-AL*V(k-INTER_Y)-AP*V(k-1)+A*V(k)-AP*V(k+1)-AL*V(k+INTER_Y)
             ELSE!BC
                E(k)=-AL*V(k-INTER_Y)-AP*V(k-1)+A*V(k)-AL*V(k+INTER_Y)
             END IF
          END DO
       ELSE!BC
          DO j=GAP(i,1),GAP(i,2)
             k=j+(i-1)*INTER_Y
             IF(j==1)THEN!BC
                E(k)=-AL*V(k-INTER_Y)+A*V(k)-AP*V(k+1)
             ELSE IF(j<INTER_Y)THEN!BC
                E(k)=-AL*V(k-INTER_Y)-AP*V(k-1)+A*V(k)-AP*V(k+1)
             ELSE!BC
                E(k)=-AL*V(k-INTER_Y)-AP*V(k-1)+A*V(k)
             END IF
          END DO
       END IF
    END DO
  END SUBROUTINE PMV
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
!COMMUNICATION VERSION 2: 2_PtoP: POINT A POINT RANG INF AU PIVOT ET POINT A POINT RANG SUP AU PIVOT:
!BUG POUR NX=1002=NY, Err N2 ET Ninfinity SONT BONNES POUR DES TAILLES PLUS PETITES
SUBROUTINE COMM_VECT_PIVOT(status,stateinfo,rang,NBR_P,ES,it1,itN,INTER_Y,INTER_X,mode,PIVOT)
    INTEGER,DIMENSION(MPI_STATUS_SIZE),INTENT(inout)::status
    INTEGER,INTENT(INOUT)::stateinfo
    INTEGER,INTENT(IN)::rang,NBR_P,it1,itN,INTER_Y,INTER_X,mode,PIVOT
    REAL(RK_GC),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(INOUT)::ES!LE VECTEUR A COMMUNIQUER
    INTEGER::NP,NPP,T1M,T1MM,REQUEST,A,B
    NP=itN+1;NPP=itN+INTER_Y;T1M=it1-1;T1MM=it1-INTER_Y
    A=itN-INTER_Y+1;B=it1+INTER_Y-1
       IF(NBR_P<=INTER_X)THEN
          IF(rang==0)THEN
             CALL MPI_SEND(ES(A:itN),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
             CALL MPI_RECV(ES(NP:NPP),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo) 
             !SEND PARTIE NORD DE ES AU RANG+1 ET RECV PARTIE SUD DE ES DU RANG+1

          ELSE IF(rang>0.AND.rang<PIVOT)THEN
             CALL MPI_RECV(ES(T1MM:T1M),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
             CALL MPI_SEND(ES(it1:B),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
             !RECV PARTIE NORD DE ES DU RANG-1 ET SEND PARTIE SUD DE ES AU RANG-1 
             CALL MPI_SEND(ES(A:itN),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
             CALL MPI_RECV(ES(NP:NPP),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
             !SEND PARTIE NORD DE ES AU RANG+1 ET RECV PARTIE SUD DE ES DU RANG+1                

          ELSE IF(rang>=PIVOT.AND.rang<NBR_P-1)THEN
             CALL MPI_RECV(ES(NP:NPP),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
             CALL MPI_SEND(ES(A:itN),INTER_Y,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)

             CALL MPI_SEND(ES(it1:B),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
             CALL MPI_RECV(ES(T1MM:T1M),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
             
             !RECV PARTIE NORD DE ES DU RANG-1 ET SEND PARTIE SUD DE ES AU RANG-1
          ELSE!RANG=NBR_P-1
             CALL MPI_SEND(ES(it1:B),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
             CALL MPI_RECV(ES(T1MM:T1M),INTER_Y,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
          END IF
          
       ELSE
          CALL COMM_NP_SUP_PIVOT(NBR_P,rang,it1,itN,INTER_Y,ES,stateinfo,status,PIVOT)!SI Np>INTER_X
       END IF
    
  END SUBROUTINE COMM_VECT_PIVOT
  !*********** UNIQUEMENT SI Np>INTER_X: ***********************************************
  subroutine COMM_NP_SUP_PIVOT(Np,rang,it1,itN,INTER,ENTER,stateinfo,status,PIVOT)
    integer::k,wr!i,j,k,liste
    integer,intent(in)::Np,rang,it1,itN,INTER,PIVOT!INTER=INTER_Y

    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    real(RK_GC),dimension(it1-INTER:itN+INTER),intent(inout)::ENTER

    INTEGER::A,B,C,D,E,F
    A=itN-INTER+1;B=it1-INTER;C=it1-1;D=it1+INTER-1;E=itN+1;F=itN+INTER
    k=0
    wr=0
    IF(rang<PIVOT)THEN
       !** BOUCLE COMMUNICATION POINT A POINT: RANG SEND PARTIE NORD DE ES A RANG+1 ET RANG+1 RECV PARTIE NORD DE ES  DE RANG
       do wr=0,PIVOT-2
          if(wr==rang)then
             if(rang<Np-1)then
                call MPI_SEND(ENTER(A:itN),INTER,MPI_DOUBLE_PRECISION,rang+1,21,MPI_COMM_WORLD,stateinfo)
                !call MPI_SEND(ENTER(itN-INTER+1:itN),INTER,MPI_DOUBLE_PRECISION,rang+1,21,MPI_COMM_WORLD,stateinfo)
             end if
          else if(rang==wr+1)then
             call MPI_RECV(ENTER(B:C),INTER,MPI_DOUBLE_PRECISION,rang-1,21,MPI_COMM_WORLD,status,stateinfo)
             !call MPI_RECV(ENTER(it1-INTER:it1-1),INTER,MPI_DOUBLE_PRECISION,rang-1,21,MPI_COMM_WORLD,status,stateinfo)
          end if
       end do
    ELSE
       !** BOUCLE COMMUNICATION POINT A POINT: SENS INVERSE
       do wr=Np-1,PIVOT+1,-1
          if(wr==rang)then
             if(rang>0)then
                call MPI_SEND(ENTER(it1:D),INTER,MPI_DOUBLE_PRECISION,rang-1,31,MPI_COMM_WORLD,stateinfo)
                !call MPI_SEND(ENTER(it1:it1+INTER-1),INTER,MPI_DOUBLE_PRECISION,rang-1,31,MPI_COMM_WORLD,stateinfo)
             end if
          else if(rang==wr-1)then
             call MPI_RECV(ENTER(E:F),INTER,MPI_DOUBLE_PRECISION,rang+1,31,MPI_COMM_WORLD,status,stateinfo)
             !call MPI_RECV(ENTER(itN+1:itN+INTER),INTER,MPI_DOUBLE_PRECISION,rang+1,31,MPI_COMM_WORLD,status,stateinfo)
          end if
       end do
    END IF
    PRINT*,'ECHO1',RANG
    !COMMUNICATION(ECHANGE) ENTRE PIVOT ET PIVOT-1:
    IF(rang==PIVOT)THEN
       CALL MPI_SEND(ENTER(it1:D),INTER,MPI_DOUBLE_PRECISION,rang-1,31,MPI_COMM_WORLD,stateinfo)
       call MPI_RECV(ENTER(B:C),INTER,MPI_DOUBLE_PRECISION,rang-1,21,MPI_COMM_WORLD,status,stateinfo)
    ELSE IF(rang==PIVOT-1)THEN
       CALL MPI_RECV(ENTER(E:F),INTER,MPI_DOUBLE_PRECISION,rang+1,31,MPI_COMM_WORLD,status,stateinfo)
       call MPI_SEND(ENTER(A:itN),INTER,MPI_DOUBLE_PRECISION,rang+1,21,MPI_COMM_WORLD,stateinfo)
    END IF
    PRINT*,'ECHO2',RANG
    CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
    PRINT*,'ECHO3',RANG
    !INVERSION:
     IF(rang>=PIVOT)THEN
       !** BOUCLE COMMUNICATION POINT A POINT: RANG SEND PARTIE NORD DE ES A RANG+1 ET RANG+1 RECV PARTIE NORD DE ES  DE RANG
       do wr=PIVOT,Np-2
          if(wr==rang)then
             if(rang<Np-1)then
                call MPI_SEND(ENTER(A:itN),INTER,MPI_DOUBLE_PRECISION,rang+1,21,MPI_COMM_WORLD,stateinfo)
                !call MPI_SEND(ENTER(itN-INTER+1:itN),INTER,MPI_DOUBLE_PRECISION,rang+1,21,MPI_COMM_WORLD,stateinfo)
             end if
          else if(rang==wr+1)then
             call MPI_RECV(ENTER(B:C),INTER,MPI_DOUBLE_PRECISION,rang-1,21,MPI_COMM_WORLD,status,stateinfo)
             !call MPI_RECV(ENTER(it1-INTER:it1-1),INTER,MPI_DOUBLE_PRECISION,rang-1,21,MPI_COMM_WORLD,status,stateinfo)
          end if
       end do
    ELSE
       !** BOUCLE COMMUNICATION POINT A POINT: SENS INVERSE
       do wr=PIVOT-1,1,-1
          if(wr==rang)then
             if(rang>0)then
                call MPI_SEND(ENTER(it1:D),INTER,MPI_DOUBLE_PRECISION,rang-1,31,MPI_COMM_WORLD,stateinfo)
                !call MPI_SEND(ENTER(it1:it1+INTER-1),INTER,MPI_DOUBLE_PRECISION,rang-1,31,MPI_COMM_WORLD,stateinfo)
             end if
          else if(rang==wr-1)then
             call MPI_RECV(ENTER(E:F),INTER,MPI_DOUBLE_PRECISION,rang+1,31,MPI_COMM_WORLD,status,stateinfo)
             !call MPI_RECV(ENTER(itN+1:itN+INTER),INTER,MPI_DOUBLE_PRECISION,rang+1,31,MPI_COMM_WORLD,status,stateinfo)
          end if
       end do
    END IF
  end subroutine COMM_NP_SUP_PIVOT
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************
! *********************************************************************************************************************************


!PGC NON UTILISE ICI:
SUBROUTINE PGC(A,AP,AL,U,F,F_VAR,userchoice,DT,TOL,KMAX,it1,itN&
       ,INTER_Y,INTER_X,S1,S2,GAP,Np,rang,stateinfo,status,mode,r,z,p)
    REAL(RK_GC),INTENT(IN)::A,AP,AL
    INTEGER,INTENT(IN)::it1,itN,INTER_Y,INTER_X,S1,S2,Np,rang,userchoice,KMAX,mode
    INTEGER,INTENT(IN),DIMENSION(S1:S2,1:2)::GAP
    REAL(RK_GC),INTENT(IN)::DT,TOL,F_VAR
    REAL(RK_GC),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(INOUT)::U
    REAL(RK_GC),DIMENSION(it1:itN),INTENT(IN)::F
    INTEGER,DIMENSION(MPI_STATUS_SIZE),INTENT(inout)::status
    INTEGER,INTENT(INOUT)::stateinfo
    REAL(RK_GC)::beta,beta_RED,alpha,beta_old,DOTzp_RED,DOTzp,gamma
    REAL(RK_GC),DIMENSION(it1:itN),INTENT(INOUT)::r,z
    REAL(RK_GC),DIMENSION(it1-INTER_Y:itN+INTER_Y),INTENT(INOUT)::p
    REAL(RK_GC),DIMENSION(it1:itN)::mp
    INTEGER::k
    !ALLOCATE(r(it1:itN),z(it1:itN),p(it1-INTER_Y:itN+INTER_Y))
    !INITIALISATION:
    k=0
    r=0.0d0;z=0.0d0
    !UN PRODUIT MATRICE/VECTEUR NECESSITE UNE COMM DE PARTIE DU VECTEUR ENTRE PROCS:
    IF(Np>1)THEN
       CALL COMM_VECT(status,stateinfo,rang,Np,U,it1,itN,INTER_Y,INTER_X,mode)
    END IF
    CALL PMV(A,AP,AL,INTER_Y,INTER_X,it1,itN,r,U,S1,S2,GAP)
    !ON NE SOUHAITE PAS RECALCULER F3, NI STOCKER b=U+DT*F
    SELECT CASE(userchoice)
    CASE(3)!POUR F=F3, DANS CE CAS ON RAJOUTE LA VARIABLE EN TEMPS POUR LA FORCE F_VAR=cos(Pi*t/2)
       r=-r+U(it1:itN)+DT*F_VAR*F!.EQV. r=b-A.U avec b=U(it1:itN)+DT*F_VAR*F
    CASE DEFAULT!POUR F=F1 OU F=F2
       r=-r+U(it1:itN)+DT*F!.EQV. r=b-A.U avec b=U(it1:itN)+DT*F
    END SELECT
    
    mp=r/A!LE PRECONDITIONNEUR DIAGONALE

    beta_RED=DOT_PRODUCT(r,mp)!NORME2 AU CARRE DU RESIDU (PAR MORCEAU)
    CALL MPI_ALLREDUCE(beta_RED,beta,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)!reduction
    

    p(it1-INTER_Y:it1-1)=0.0d0;p(it1:itN)=mp(it1:itN);p(itN+1:itN+INTER_Y)=0.0d0!p la direction de descente du GC
    
    DO WHILE(SQRT(beta)>TOL.AND.(k<KMAX))!BOUCLE DE CONVERGENCE
       !z=Ap, OBLIGATION DE FAIRE UNE COMM SUR p:
       IF(Np>1)THEN
          CALL COMM_VECT(status,stateinfo,rang,Np,p,it1,itN,INTER_Y,INTER_X,mode)
       END IF
       CALL PMV(A,AP,AL,INTER_Y,INTER_X,it1,itN,z,p,S1,S2,GAP)!PRODUIT MATRICE VECTEUR

       !REDUCTION SUR DOT_PRODUCT(z,p):
       DOTzp_RED=DOT_PRODUCT(z(it1:itN),p(it1:itN))
       CALL MPI_ALLREDUCE(DOTzp_RED,DOTzp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
       
       !CALCUL DE alpha le pas de descente:
       alpha=beta/DOTzp

       !ACTUALISATION DE U ET de r:
       U(it1:itN)=U(it1:itN)+alpha*p(it1:itN)
       r=r-alpha*z

       !SAVE beta dans beta_old:
       beta_old=beta
       
       !ACTUALISATION DE mp:
       mp=r/A
       
       !ACTUALISATION DE BETA ET REDUCTION:
       beta_RED=DOT_PRODUCT(r,mp)
       CALL MPI_ALLREDUCE(beta_RED,beta,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,stateinfo)
       
       !CALCUL de gamma:
       gamma=beta/beta_old
       
       !ACTUALISATION DE LA DIRECTION DE DESCENTE p:
       p(it1:itN)=mp(it1:itN)+gamma*p(it1:itN)

       !ACTUALISATION DE k:
       k=k+1
              
    END DO
    PRINT*,'k',k,'SQRT(beta)',SQRT(beta),'TOL',TOL,'RANG',rang
    
  END SUBROUTINE PGC






end module ModGC
