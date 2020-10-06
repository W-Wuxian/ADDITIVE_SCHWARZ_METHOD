module ModData
  implicit none

  !Nx-2 noeuds interieurs suivant x
  !Ny-2 noeuds interieurs suivant y
  INTEGER,parameter::RK_DATA=8
  !ModDATA CONTIENT LA SUBROUTINE POUR LIRE LE FICHIER DATA REMPLIT PAR L'UTILISATEUR
contains

  subroutine INI_PARA(rang,Lx,Ly,D,Nx,Ny,PARAMETRES_USER,sysmove,userchoice,DT,ITER_TMAX)
    integer,intent(in)::rang
    !Lx, Ly longueur suivant x et y
    !D coeff de diffusion
    !Nx Ny nombre de noeuds dans la direction x et y
    !ch nom du fichier à lire pour obtenir les para
    real(RK_DATA),intent(out)::Lx,Ly,D,DT
    integer,intent(out)::Nx,Ny,sysmove,userchoice,ITER_TMAX
    character(len=4),intent(in)::PARAMETRES_USER

    open(10+rang,file=PARAMETRES_USER)
    read(10+rang,*)
    read(10+rang,*)Lx,Ly,D,Nx,Ny!LX,Ly,Nx,Ny:PARAMETRES GEOMETRIQUES ET D: COEFF DIFFUSION 
    read(10+rang,*)
    read(10+rang,*)sysmove!=1 => VISUALISATION SOLUTION; =0 => PAS DE VISUALISATION, JUSTE ECRITURE SOLUTION
    read(10+rang,*)
    read(10+rang,*)
    read(10+rang,*)
    read(10+rang,*)
    read(10+rang,*)
    read(10+rang,*)userchoice!LE CAS SOURCE:1=F1;2=F2;3=F3
    read(10+rang,*)
    read(10+rang,*)DT,ITER_TMAX!LE PAS DE TEMPS ET L'ITERATION MAX EN TEMPS
    close(10+rang)
    

  end subroutine INI_PARA

end module ModData
