 !-----------------------------------------
! CE MODULE CONTIENT LE GRADIENT CONJUGUÉ
! EN PLEIN ET EN CREUX (FORMAT COORDONNES)
!-----------------------------------------

Module mod_gradconj
!---modules----------------
  use mod_parametres
  use mod_fonctions
!--------------------------
  Implicit None


Contains

!---GRADIENT CONJUGUÉ POUR UNE MATRICE A PLEINE---------
  Subroutine gradconj(A,b,x0,x,beta,k,n)
    !---variables------------------------------------
    Integer, Intent(In) :: n                     !taille vecteur solution (=taille de la matrice carrée)
    Real(PR), Dimension(:,:), Intent(In) :: A    !matrice à inverser
    Real(PR), Dimension(:), Intent(In) :: b,x0   ! b second membre, x0 CI
    Real(PR), Dimension(:), Intent(Out) :: x     !solution finale
    Real(PR), Intent(Out) :: beta                !norme du résidu
    Real(PR), Dimension(:), Allocatable :: p,z,r1,r2
    Real(PR) :: alpha,gamma,eps
    Integer, Intent(Out) :: k                   !nb d'itération pr convergence
    Integer :: kmax
    !-------------------------------------------------
    eps=0.01_PR
    kmax=150
    Allocate(p(n),z(n),r1(n),r2(n))
    r1=b-matmul(A,x0)
    p=r1
    beta=sqrt(sum(r1*r1))
    k=0
    Do While (beta>eps .and. k<=kmax)
       z=matmul(A,p)
       alpha=dot_product(r1,r1)/dot_product(z,p)
       x=x+alpha*p
       r2=r1-alpha*z
       gamma=dot_product(r2,r2)/dot_product(r1,r1)
       p=r2+gamma*p
       beta=sqrt(sum(r1*r1))
       k=k+1
       r1=r2
    End Do
    If (k>kmax) then
       Print*,"Tolérance non atteinte:",beta
    End If
    Deallocate(p,z,r1,r2)
  End Subroutine gradconj
!----------------------------------------------------

!---GRADIENT CONJUGUÉ POUR MATRICE CREUSE AU FORMAT COORDONNÉES-----
!---A PARALLELISER : PRODUIT MATRICE/VECTEUR + PRODUIT SCALAIRE
 Subroutine gradconj_SPARSE(AA,IA,JA,F,x0,x,b,k,n)
    Implicit None
    !---variables ---------------------------------------
    Integer, Intent(In) :: n   !taille vecteur solution
    Real(PR), Dimension(:), Intent(In) :: AA
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(:), Intent(In) :: f,x0  !f=vect source , x0=CI
    Real(PR), Dimension(:), Intent(InOut) :: x  !solution finale
    Real(PR), Intent(Out) :: b  !norme du résidu
    Real(PR), Dimension(:), Allocatable :: p,z,r1,r2
    Real(PR) :: a,g,eps
    Integer, Intent(Out) :: k !nb d'itérations pour convergence
    Integer :: kmax   !nb d'itérations max
    !-----------------------------------------------------
    eps=0.01_PR
    kmax=1500
    Allocate(p(n),z(n),r1(n),r2(n))
    !---paralléliser------------------------
    r1=f-MatVecSPARSE(AA,IA,JA,x0)
    !-------------------------------------
    p=r1
    b=sqrt(sum(r1*r1))
    k=0
    Do While (b>eps .and. k<=kmax)
       !---paralléliser------------------
       z=MatVecSPARSE(AA,IA,JA,p)
       !---------------------------------
       a=dot_product(r1,r1)/dot_product(z,p)
       x=x+a*p
       r2=r1-a*z
       g=dot_product(r2,r2)/dot_product(r1,r1)
       p=r2+g*p
       b=sqrt(sum(r1*r1))
       k=k+1
       r1=r2
    End Do
    If (k>kmax) then
       Print*,"Tolérance non atteinte:",b
    End If
    Deallocate(p,z,r1,r2)
  End Subroutine gradconj_SPARSE
!--------------------------------------------------------------

!---FONCTION PRODUIT MATxVEC EN SPARSE-------------------------
 Function MatVecSPARSE(AA,IA,JA,x) Result(y)
    !Produit MatriceSPARSE.vecteur plein (x) retourne vecteur plein(y)
    !---variables------------------------------------------
    Real(PR), Dimension(:), Intent(In) :: AA,x
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(Size(x)) :: y
    Integer :: i,n,nnz
    !--------------------------------------------------------
    n=Size(x,1)
    nnz=Size(AA,1)
    y(1:n)=0._PR
    Do i=1,nnz
          if (IA(i)==0) then
            print*,"element IA nul pour i =",i
          end if
          y(IA(i))=y(IA(i))+AA(i)*x(JA(i))
    End Do
  End Function MatVecSPARSE
!---------------------------------------------------------------







!GC COO FORMAT BUILD  FOR global numerotation over x
Subroutine MOD_gradconj_SPARSE(AA,IA,JA,f,x0,x,b,k,it1,itN)
    Implicit None
    !---variables ---------------------------------------
    Integer, Intent(In) :: it1,itN   !taille vecteur solution
    Real(PR), Dimension(:), Intent(In) :: AA
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(it1:itN), Intent(In) :: f,x0  !f=vect source , x0=CI
    Real(PR), Dimension(it1:itN), Intent(InOut) :: x  !solution finale
    Real(PR), Intent(Out) :: b  !norme du résidu
    Real(PR), Dimension(:), Allocatable :: p,z,r1,r2
    Real(PR) :: a,g,eps
    Integer, Intent(Out) :: k !nb d'itérations pour convergence
    Integer :: kmax   !nb d'itérations max
    !-----------------------------------------------------
    eps=0.0000000001_PR
    kmax=2*Na_l!1500
    Allocate(p(it1:itN),z(it1:itN),r1(it1:itN),r2(it1:itN))
    !---paralléliser------------------------
    r1=f-MOD_MatVecSPARSE(AA,IA,JA,x0,it1,itN)
    !-------------------------------------
    p=r1
    b=sqrt(sum(r1*r1))
    k=0
    Do While (b>eps .and. k<=kmax)
       !---paralléliser------------------
       z=MOD_MatVecSPARSE(AA,IA,JA,p,it1,itN)
       !---------------------------------
       a=dot_product(r1,r1)/dot_product(z,p)
       x=x+a*p
       r2=r1-a*z
       g=dot_product(r2,r2)/dot_product(r1,r1)
       p=r2+g*p
       b=sqrt(sum(r1*r1))
       k=k+1
       r1=r2
       !PRINT*,'k,b',k,b
    End Do
    If (k>kmax) then
       Print*,"Tolérance non atteinte:",b
    End If
    Deallocate(p,z,r1,r2)
  End Subroutine MOD_gradconj_SPARSE
!--------------------------------------------------------------



!---FONCTION PRODUIT MATxVEC EN SPARSE-------------------------
 Function MOD_MatVecSPARSE(AA,IA,JA,x,it1,itN) Result(y)
    !Produit MatriceSPARSE.vecteur plein (x) retourne vecteur plein(y)
    !---variables------------------------------------------
   INTEGER,INTENT(IN)::it1,itN
    Real(PR), Dimension(:), Intent(In) :: AA
    Real(PR), Dimension(it1:itN), Intent(In) ::x
    Integer, Dimension(:), Intent(In) :: IA,JA
    
    Real(PR), Dimension(it1:itN) :: y
    Integer :: i,n,nnz,val
    !--------------------------------------------------------
    n=Size(x,1)
    !PRINT*,'***',n,itN-it1+1
    nnz=Size(AA,1)
    y(it1:itN)=0._PR
    val=(it1-1)
    Do i=1,nnz
          if (IA(i)==0) then
            print*,"element IA nul pour i =",i
          end if
          !PRINT*,i,IA(i),val+IA(i),JA(i),val+JA(i),x(val+JA(i))
          y(val+IA(i))=y(val+IA(i))+AA(i)*x(val+JA(i))
          
    End Do
  End Function MOD_MatVecSPARSE
!---------------------------------------------------------------




Subroutine MOD_jacobi_SPARSE(AA,IA,JA,F,x0,x,resnorm,k)
 
  !-----variables---------------------
  
  Real(PR), Dimension(nnz_l), Intent(In) :: AA
  Integer, Dimension(nnz_l), Intent(In) :: IA,JA
  Real(PR), Dimension(it1:itN), Intent(InOut) :: F,x0  !f=vect source , x0=CI
  Real(PR), Dimension(it1:itN), Intent(InOut) :: x  !solution finale
  Real(PR), Intent(Out) :: resnorm  !norme du résidu
  Real(PR), Dimension(:), Allocatable :: r
  Real(PR) :: eps,somme
  Integer, Intent(InOut) :: k !nb d'itérations pour convergence
  Integer :: kmax,pp,val,ii,borne   !nb d'itérations max
  ALLOCATE(r(it1:itN))
  
  borne=nnz_l+1
  x=x0

  eps = 0.0000000001_PR
  kmax = 2*Na_l*Na_l
  r = F-MOD_MatVecSPARSE(AA,IA,JA,x0,it1,itN) !résidu initial
  resnorm = sqrt(sum(r*r)) !norme 2 du résidu
  k = 0 !initialisation itération
  boucle_resolution: Do While ((resnorm>eps) .and. (k<kmax))

    k=k+1
    !on calcule xi^k+1=(1/aii)*(bi-sum(j=1,n;j/=i)aij*xj^k)
    x0=0._PR
    pp=1;val=it1-1
    
    boucle_xi: Do ii=1,Na_l
       !calcul de la somme des aij*xj^k pour j=1,n;j/=i
       somme = 0.0_PR;
!!$       if(rang==0)then
!!$          PRINT*,'ii,pp,k,nnz_l,ubound IA',ii,pp,k,nnz_l,ubound(IA,1),nnz_l+1
!!$       end if
      
          DO WHILE((pp<borne) .AND. (ii==IA(pp)))
            
                IF (IA(pp) /= JA(pp)) THEN

                   somme= somme + x(val+JA(pp)) * AA(pp)
                END IF

!!$             if(rang==0)then
!!$                PRINT*,'ii,pp,IA(pp),JA(pp),AA(pp),x(JA(pp))',ii,pp,IA(pp),JA(pp),AA(pp),x(val+JA(pp))
!!$             end if
                pp=pp+1
                !PRINT*,'pp',pp
                if(pp==borne)then !because the do while check both pp<borne and ii==IA(pp),maybe due to makefile
                   ! he should break when pp>=borne and don't check ii==IA(pp)
                   GOTO 7777
                end if
          END DO
      
      
       !calcul du nouveau xi
7777      x0(val+ii)=(1._PR/alpha)*(F(val+ii)-somme)
          x0(val+ii)=(1._PR/alpha)*(F(val+ii)-somme)!;print*,'t',rang
    End Do boucle_xi
    x=x0
    !on calcule le nouveau résidu
    r = F - MOD_MatVecSPARSE(AA,IA,JA,x,it1,itN)
    ! calcul de la norme 2 du résidu
    resnorm = sqrt(sum(r*r))
    !print*,resnorm,k

  End Do boucle_resolution
  If (k>kmax) then
     Print*,"Tolérance non atteinte:",resnorm
  End If
  DEALLOCATE(r)
  
End Subroutine MOD_jacobi_SPARSE


Subroutine BICGSTAB(it1,itN,AA,IA,JA,x,b,eps)
  Implicit None
  Integer, intent(IN) :: it1, itN
  Real(PR), Dimension(:), intent(IN) :: AA
  Integer,  Dimension(:), intent(IN) :: IA, JA
  Real(PR), Dimension(:), intent(INOUT) :: x
  Real(PR), Dimension(:), intent(IN) :: b
  Real(PR), intent(IN) :: eps
  Real(PR), Dimension(it1:itN) :: r0, v, p ,r, h, s, t
  Real(PR) :: rho_i, rho_im1, alpha_CG, w, beta_cg
  Real(PR) :: omega
  Real(PR) :: beta
  Integer  :: k, kmax

  kmax  = 2*Na_l*Na_l
  omega = 1.0_PR

  r0 = b -MOD_MatVecSPARSE(AA, IA, JA, x, it1, itN)
  r  = r0
  rho_im1 = 1.0_PR; alpha_CG = 1.0_PR; w = 1.0_PR
  v = 0.0_PR; p = 0.0_PR

  k = 0
  beta = SQRT(DOT_PRODUCT(r0,r0))

  DO WHILE(beta>eps.AND.k<kmax)
     rho_i = DOT_PRODUCT(r0,r)
     beta_CG = (rho_i/rho_im1)*(alpha_CG/w)
     p = r + beta_CG*(p-v*w)
     v = MOD_MatVecSPARSE(AA,IA,JA,p,it1,itN)
     alpha_CG = rho_i/DOT_PRODUCT(r0,v)
     h = x + p*alpha_CG
     s = r - v*alpha_CG
     t = MOD_MatVecSPARSE(AA,IA,JA,s,it1,itN)
     w = DOT_PRODUCT(t,s)/DOT_PRODUCT(t,t)
     x = h + s*w
     r = s - t*w
     beta = DOT_PRODUCT(r,r)
     rho_im1 = rho_i
     k  = k+1
  END DO
  !PRINT*,'kmax,k',kmax,k
 
End Subroutine BICGSTAB




End Module mod_gradconj
