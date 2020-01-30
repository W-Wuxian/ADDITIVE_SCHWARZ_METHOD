module mod_gmres
  USE mod_parametres
  USE mod_fonctions
Implicit None


Contains

  !---FONCTION PRODUIT MATxVEC EN SPARSE-------------------------
   Function MatVecSPARSE(AA,IA,JA,x) Result(y)
     Implicit None
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
              !print*,"element IA nul pour i =",i
            end if
            y(IA(i))=y(IA(i))+AA(i)*x(JA(i))

      End Do
    End Function MatVecSPARSE
  !---------------------------------------------------------------
  !---FONCTION PRODUIT MATxVEC EN SPARSE-------------------------
  Function COL_MatVecSPARSE(AA,IA,JA,x,dim) Result(y)
    Implicit None
    !Produit MatriceSPARSE.vecteur plein (x) retourne vecteur plein(y)
    !---variables------------------------------------------
    Integer, Intent(IN) :: dim
    Real(PR), Dimension(:), Intent(In) :: AA
    Real(PR), Dimension(1:dim,1), Intent(In) :: x
    Integer, Dimension(:), Intent(In) :: IA,JA
    Real(PR), Dimension(dim) :: y
    Integer :: i,nnz
    !--------------------------------------------------------
    !n=Size(x,1)
    nnz=Size(AA,1)
    y(1:dim)=0._PR
    Do i=1,nnz
       if (IA(i)==0) then
          !print*,"element IA nul pour i =",i
       end if
       y(IA(i))=y(IA(i))+AA(i)*x(JA(i),1)

    End Do
  End Function COL_MatVecSPARSE
  !---------------------------------------------------------------

  Subroutine GMres_COO(AA,IA,JA,f,x0,x,b,k,dim,it1,itN,DIMKRYLOV)
    Implicit None
    !-----Variables: dim = Na_l
    Integer, Intent(IN)                         :: dim, it1, itN, DIMKRYLOV
    Real(PR), Dimension(:), Intent(In)          :: AA
    Integer, Dimension(:), Intent(In)           :: IA,JA
    Real(PR), Dimension(it1:itN), Intent(In)    :: f,x0  !f=vect source , x0=CI
    Real(PR), Dimension(it1:itN), Intent(InOut) :: x  !solution finale
    Real(PR), Intent(Out)                       :: b  !norme du résidu
    Integer, Intent(Out)                        :: k  !nb d'itérations pour convergence

    Integer  :: m      ! m = DIMKRYLOV
    Integer  :: kmax   !nb d'itérations max
    Real(PR) ::beta    !norm2 de r
    Real(PR) :: eps    !tol solver
    ! Rescale: Bufferx(1:dim) = x(it1:itN) so natural dot_product
    Real(PR), Dimension(:), Allocatable :: Bufferx0, Bufferx, Bufferf
    !---- Vectors & Matrix used for GMres:
    Real(PR),dimension(1:dim,1:DIMKRYLOV+1)::Vm
    !Real(PR),dimension(1:m+1,1:m)::Hm_
    Real(PR),dimension(:,:),allocatable::Hm_
    Real(PR),dimension(1:DIMKRYLOV+1,1:DIMKRYLOV)::Rm_
    Real(PR),dimension(1:DIMKRYLOV+1,1:DIMKRYLOV+1)::Qm_
    Real(PR),dimension(1:dim,1:DIMKRYLOV)::FF !Vm+1.Hm_
    Real(PR),dimension(1:DIMKRYLOV+1)::gm_
    Real(PR),dimension(1:DIMKRYLOV)::y !sol de Rm.y=gm
    Real(PR),dimension(1:dim)::Vmy
    Real(PR),dimension(1:dim)::AX
    Real(PR),dimension(1:dim)::z,p,r


    Allocate(Bufferx0(1:itN-it1+1))
    Bufferx0(1:dim) = x0(it1:itN)
    Allocate(Bufferx(1:itN-it1+1))
    Bufferx(1:dim) = x(it1:itN)
    Allocate(Bufferf(1:itN-it1+1))
    Bufferf(1:dim) = f(it1:itN)

    m = DIMKRYLOV
    beta=0.0_PR
    k=0

    !**INITIALISATION DE ro=b-AXo:MatVecSPARSE(AA,IA,JA,x)
    r=Bufferf-MatVecSPARSE(AA,IA,JA,Bufferx0)
    beta = DOT_PRODUCT(r,r)
    beta = sqrt(beta)

    eps = 0.00001_PR
    kmax = dim**2

    boucle_eps_kmax:DO WHILE(beta>eps .AND. k<kmax)

       !Construction DE LA BASE Vm et de Hm_
       ALLOCATE(Hm_(1:m+1,1:m))
       CALL arnoldi_reortho(dim,m,AA,IA,JA,r,beta,Vm,Hm_)
       !Decomposion QR de Hm_ pour avoir Qm_ et Rm_
       Rm_ = Hm_
       DEALLOCATE(Hm_)
       CALL givens_QR_opt2(dim,m,Rm_,Qm_)
       !CALL givens_QR_opt(dim,m,Hm_,Qm_,Rm_)
       !Setup gm_ = transpose(Qm_). beta*e1
       gm_LOOP:DO i = 1,m+1
          gm_(i) = beta * Qm_(1,i)
       END DO gm_LOOP
       !Resolution DE Rm.y=gm où Rm=Rm_(1:m,:)et gm= gm_(1:m)
       y(m)=gm_(m)/Rm_(m,m)
       moindrecarre:do j=m-1,1,-1
          y(j)=gm_(j)
          do c=j+1,m
             y(j)=y(j)-Rm_(j,c)*y(c)
          end do
          y(j)=y(j)/Rm_(j,j)
       end do moindrecarre
       !Calcul de Vmy
       Vmy=MATMUL(Vm(1:dim,1:m),y)
       !CALCUL de X=X+Vm.y
       Bufferx=Bufferx+Vmy
       !Actualisation r:
       r=bufferf-MatVecSPARSE(AA,IA,JA,Bufferx)
       beta=sqrt(DOT_PRODUCT(r,r))
       k=k+1
      ! PRINT*,'APRES k=k+1',k,beta
       if(rang==0)then
          print*,'gmres rang,k,res',rang,k,beta
       end if
    END DO boucle_eps_kmax

    b = beta
    x(it1:itN) = Bufferx(1:dim)
    Deallocate(Bufferx0,Bufferx,Bufferf)

  End Subroutine GMres_COO


  !----------
  !$$$$$$$$$$$$$$$$ARNOLDI MGS REORTHO$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  subroutine arnoldi_reortho(dim,m,AA,IA,JA,ro,beta,Vm,Hm_)!v,Vm,Hm_)
    Implicit None
    !la notation Hm_  pour la barre up du Hm correspondant
    !à la notation theorique
    !m nbr de colonne de Vm ; Vm appartient à M_n,m(R), n =dim ici
    integer,intent(in)::dim,m
    !on utilise m pour fixer la dimension de Vm,
    !base constituée des (vi)i
    Real(PR),dimension(1:dim),intent(in)::ro
    Real(PR),intent(in)::beta
    !ro=b-Axo
    !beta =||ro|| au carré
    !Real(PR),dimension(1:dim),intent(inout)::v
    !v vecteur constituant les colonnes de Vm

    Real(PR), Dimension(:), Intent(In) :: AA
    Integer, Dimension(:), Intent(In) :: IA,JA
    !Real(PR),dimension(1:dim,1:dim),intent(in)::A


    Real(PR),dimension(1:dim,1:m+1),intent(out)::Vm
    Real(PR),dimension(1:m+1,1:m),intent(out)::Hm_


    Real(PR),dimension(1:dim)::z,Av
    !z vecteur de stockage pour calcul de v
    !Av pour stocker A.v_j ,v_j la  colonne de Vm
    Real(PR)::pscal,nz,tmp !nz pour calculer la norme de z
    integer::k,i,j,c,d

    !**INITIALISATION DE v:
    Vm=0.0_PR
    Hm_=0.0_PR
    !v=ro/beta !v1
    do i=1,dim
       Vm(i,1)=ro(i)/beta
    end do


    !**CALCUL DE z:
    j=0
    do while(j<m)
       j=j+1
       !$$CALCUL DE A.v:
       Av=0.0_PR

       Av = Av + COL_MatVecSPARSE(AA,IA,JA,Vm(1:dim,j),dim)

       z=Av
       do i=1,j
          !$$CALCUL DE SUM(<Avj|vi>.vi)
          pscal=0.0_PR
          pscal=DOT_PRODUCT(Av,Vm(1:dim,i))

          Hm_(i,j)=pscal

          z=z-Hm_(i,j)*Vm(:,i)

       end do!fin do i
       !$$$$$REORTHO
       do i=1,j!j
          tmp=0.0_PR
          tmp=DOT_PRODUCT(z,Vm(:,i))
          z=z-tmp*Vm(:,i)
          Hm_(i,j)=Hm_(i,j)+tmp
       end do

       !$$CALCUL DE ||z||:
       nz=0.0_PR
       nz=DOT_PRODUCT(z,z)

       nz=sqrt(nz)
       Hm_(j+1,j)=nz


       !$$CONDITION ARRET HEUREUX:
       if(abs(Hm_(j+1,j))<0.00000001_PR)then

          j=m+1!stop

       else !end if
          !$$FIN CONDITION ARRET HEUREUX

       !$$ACTUALISATION/CALCUL DE Vm(:,j+1); v_(j+1) la j+1 ieme cows de Vm
          Vm(:,j+1)=z !/Hm_(j+1,j)
          Vm(:,j+1)=Vm(:,j+1)/nz
       end if

!       PRINT*,'MGS_REORTHO',j
    end do

  end subroutine arnoldi_reortho



  subroutine givens_QR_opt2(dim,m,Rm_,Qm_)!Hm_,Qm_,Rm_)
      Implicit None
      integer,intent(in)::dim,m
      !real,dimension(1:m+1,1:m),intent(in)::Hm_
      Real(PR),dimension(1:m+1,1:m),intent(inout)::Rm_
      Real(PR),dimension(1:m+1,1:m+1),intent(out)::Qm_

      Real(PR)::c,s
      Real(PR)::coef1,coef2,coef_a,coef_b
      Real(PR)::TQ1,TQ2,TQ3,TQ4
      Real(PR)::Qg,Qd
      integer::i,j,k,l
      !Rm_=Hm_

      Qm_=0.0_PR
      do i=1,m+1
         Qm_(i,i)=1.0_PR!?stock_T_Qm_(i,i)=1.
      end do

      do j=1,m!-1!m
         do l=j+1,j+1,-1!m+1,j+1,-1  !DE j+1 à j+1 car Hm_ forme de Hessenberg



            coef1=Rm_(l,j);coef2=Rm_(l-1,j)
            c=coef1*coef1+coef2*coef2;c=sqrt(c)
            s=c
            c=coef2/c;s=-coef1/s

            TQ1=c  !T_Qm_(l,l)=c
            TQ2=c  !T_Qm_(l-1,l-1)=c
            TQ3=-s !T_Qm_(l-1,l)=-s
            TQ4=s  !T_Qm_(l,l-1)=s

            do k=j,m
               coef_a=Rm_(l-1,k);coef_b=Rm_(l,k)
               Rm_(l-1,k)=c*coef_a-s*coef_b
               Rm_(l,k)=s*coef_a+c*coef_b
            end do
            Rm_(l,j)=0.0_PR

            TQ3=-TQ3
            TQ4=-TQ4
            if(j==1 .AND. l==j+1)then
               Qm_(l,l)=TQ1;Qm_(l-1,l-1)=TQ2;Qm_(l-1,l)=TQ3;Qm_(l,l-1)=TQ4
            else
               do i=1,m+1
                  Qg=Qm_(i,l-1);Qd=Qm_(i,l)
                  Qm_(i,l-1)=Qg*TQ2+Qd*TQ4
                  Qm_(i,l)=Qg*TQ3+Qd*TQ1
               end do
            end if

         end do !fin do l
 !        PRINT*,'GIVENS',j
      end do !fin do j

    end subroutine givens_QR_opt2


end  module mod_gmres
