!-----------------------------------------------
! CE MODULE CONTIENT LES PARAMETRES DU PROBLEME
!-----------------------------------------------

Module mod_parametres
  USE mpi
  Implicit None

!---Précision des calculs -------------------
Integer, Parameter :: PR=8                       !simple précision
!--------------------------------------------

!---Parametres géométriques-----------------
Real(PR), Parameter :: Lx=1._PR  !Longueur de la barre
Real(PR), Parameter :: Ly=1._PR  !Largeur de la barre
Real(PR), Parameter :: D=1._PR   !Coefficient de diffusion
!--------------------------------------------



!PARAMETRE DE VISUALISATION:
INTEGER :: sysmove

!---Parametres numériques-------------------
!g means global Domain GLobal matrix WITHOUT DD ADDITIVE
!l means local Domain local matrix
Integer             :: nnz_g
Integer, Parameter  :: Nx_g=400              !lig nes de A (discrétisation en x)
Integer, Parameter  :: Ny_g=400           !colonnes de A (discrétisation en y)
Integer, Parameter  :: Na_g=Nx_g*Ny_g        !taille de la matrice A
Integer, Parameter  :: Nt=30                !discrétisation du temps
Real(PR), Parameter :: dx=1._PR/(Real(Nx_g)+1)     ! pas en x
Real(PR), Parameter :: dy=1._PR/(Real(Ny_g)+1)     ! pas en y
Real(PR), Parameter :: Tf=2.0_PR                 !temps final de la simulation
Real(PR), Parameter :: dt=Tf/(Real(Nt))        !pas de temps
Real(PR), Parameter :: pi=4._PR*atan(1._PR)
Real(PR), Parameter :: alpha =1._PR+(2._PR*D*dt/(dx**2._PR))+(2._PR*D*dt/(dy**2._PR)) !
Real(PR), Parameter :: beta = (-D*dt)/(dx**2._PR)      !AL                               ! CF coefficients matrice
Real(PR), Parameter :: gamma =(-D*dt)/(dy**2._PR)   !AP  
                               !
!-------------------------------------------


!---PARAMETRES MPI--------------------------
INTEGER   ::rang,Pivot
INTEGER   ::Np
INTEGER   ::stateinfo
INTEGER,DIMENSION(MPI_STATUS_SIZE)::status
CHARACTER(LEN=3) ::rank
INTEGER          ::TAG !FOR TAG FILE
!------------------------------------------

!---DEFINE INTERVALLE SUIVANT Y PAR PROCS------
INTEGER   ::S1_old, S2_old, it1_old, itN_old !avant call charge_overlap
INTEGER   ::S1
INTEGER   ::S2  !=> X_i^(rang) \in [|S1;S2|]
INTEGER   ::it1
INTEGER   ::itN !=> P(OMEGA^(rang)) \in [|it1;itN|]
!INTEGER   ::overlapd
!INTEGER   ::overlapg
INTEGER,PARAMETER::overlap=1
INTEGER   ::Na_l  !NBR rows or cols in local matrix
                  !na_loc == (S2-S1+1)*Ny
INTEGER   ::Nx_l  !Nx local will be set to S2-S1+1
INTEGER   ::Ny_l
INTEGER   ::nnz_l !nbr de non zero in local matrix
INTEGER   ::crtl_nnz_l !control de nnz
!since a A_l matrix local is made by block D AND C such:
!avec s=alpha, x=gamma et v=beta:

![D][C][0][0]     !    |s x 0 0|   !    |v 0 0 0| ! A_l.rows=Nx_l*Ny_l=A_l.cols
![C][D][C][0]     ![D]=|x s x 0|   ![C]=|0 v 0 0| ! D.rows=D.cols=Ny_l
![0][C][D][C]     !    |0 x s x|   !    |0 0 v 0| ! C.cols=C.rows=Ny_l
![0][0][C][D]     !    |0 0 x s|   !    |0 0 0 v| !We got (Nx_l) [D] block & (Nx_l-1) [C] upper block

!& (Nx_l-1) [C] lower block
!SUCH THAT crtl_nnz_l=(Nx_l*Ny_l) + 2 * (Nx_l) * (Ny_l - 1) + 2 * (Nx_l - 1) * (Ny_l)
INTEGER   ::D_rows !will be set to Ny_l
INTEGER   ::D_nnz  !will be set to D_rows+2*(D_rows-1)
INTEGER   ::C_rows !will be set to Ny_l
INTEGER   ::C_nnz  !will be set to C_rows
! WE DEFINE a control nnz_l parameter over L1, L2, L3 such that:
!|[D][C][0][0]| = [L1]
INTEGER   ::crtl_L1_nnz !will be set to D_nnz+C_nnz
!|[C][D][C][0]| = [L2] 
!|[0][C][D][C]|
INTEGER   ::crtl_L2_nnz !will be set to (Nx_l-2)*(D_nnz+2*C_nnz)
!|[0][0][C][D]| = [L3]
INTEGER   ::crtl_L3_nnz !will be set to D_nnz+C_nnz
!SUCH THAT THE RELATION (*) NEED TO BE .TRUE.:
!(*)crtl_L1_nnz+crtl_L2_nnz+crtl_L3_nnz = crtl_nnz_l = nnz_l
INTEGER   ::sum_crtl_L_nnz !sum of crtl_Li_nnz
!DANS LE CAS OU Nx_l == 1:
INTEGER   :: crtl_L1_nnz_Nxl_EQV_1 ! will be set to D_nnz

!---variables---------------------------------!
Real(PR), Dimension(:), Allocatable :: X
Real(PR), Dimension(:), Allocatable :: Y !Ny_g+2
Real(PR), Dimension(:), Allocatable :: T !Nt+2
Real(PR), Dimension(:), Allocatable :: AA
Integer,  Dimension(:), Allocatable :: IA,JA
Real(PR), Dimension(:), Allocatable :: U,Uo,F  !Na_l
Real(Pr), Dimension(:), Allocatable :: UG, UD  ! CL FROM BORDER IMMERGEE need to be set at zero
Real(PR) :: res,t1,t2              !résidu du grandconj
Integer  :: k,i,j,it,CT

INTEGER ::LBX  !init in mod_fonctions, para_param
INTEGER ::UBX 



!VARIABLES POUR LA BOUCLE SCHWARZ, CONVERGENCE DD:
REAL(PR) :: Norm1, Norm2, C1, err
REAL(PR), PARAMETER :: err_LIM = 0.0000000001_PR
INTEGER  :: step

End Module mod_parametres
