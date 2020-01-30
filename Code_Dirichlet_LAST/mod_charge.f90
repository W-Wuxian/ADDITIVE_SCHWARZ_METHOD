module mod_charge
  use mod_parametres
  Implicit None

CONTAINS
  !**************************
  !^ y Ny
  !|
  !|
  !|
  !1---------> x Nx
  !*************************
  !SE BASE SUR LA CHARGE EN X POUR CALCULER LA CHARGE TOTALE(GLOBALE)
  !INDUIT POUR MOD(DIM,Np)/=0 UNE DIFFERENCE DE CHARGE = A Ny ENTRE 2 PROCS
  !UTILISABLE POUR Np<=INTER_X
  Subroutine MODE1(rang, Np, S1, S2, it1, itN)
    Integer, Intent(IN)  :: rang, Np
    Integer, Intent(OUT) :: S1, S2,it1, itN
    !CHARGE EN X :
    CALL CHARGE_X(rang, Np, S1, S2)
    !CHARGE TOT :
    CALL CHARGE_TOT(S1, S2, it1, itN)
  End Subroutine MODE1

  !REPARTITION CLASSIQUE DE LA CHARGE EN X
  Subroutine CHARGE_X(rang, Np, S1,S2)
    Integer, Intent(IN)  :: rang, Np
    Integer, Intent(OUT) :: S1, S2
    REAL                 :: CO_RE !COEFFICIANT DE REPARTITION
    IF(Np == 1)THEN
       S1 = 1; S2 = Nx_g
    ELSE
       CO_RE=(Nx_g)/(Np)
       If(rang < mod(Nx_g,Np))Then
          S1 = rang*(CO_RE+1) + 1; S2 = (rang+1) * (CO_RE+1)
       Else
          S1 = 1 + mod(Nx_g,Np) + rang*CO_RE; S2 = S1+CO_RE-1
       End If
    END IF
  End Subroutine CHARGE_X

  !CHARGE TOTALE SUR LA NUMEROTATION GLOBALE
  Subroutine CHARGE_TOT(S1, S2, it1, itN)
    Integer, Intent(IN)  :: S1, S2
    Integer, Intent(OUT) :: it1, itN
    it1 = (S1-1)*Ny_g+1; itN = S2*Ny_g 
  End Subroutine CHARGE_TOT

  !RE-COMPUTE CHARGE WITH OVERLAP
  Subroutine CHARGE_OVERLAP(Np, rang, S1, S2, it1, itN)
    Implicit None
    Integer, Intent(IN) ::Np, rang
    Integer, Intent(INOUT) ::S1, S2, it1, itN
    !Le (-1) sur overlap pour retirer la frontiere
    !immergée qui passe en CL
   
    !WE PUT A Ny_g BECAUSE Ny_l == Ny_g and Ny_l is defined after this subroutine
    !Also Because IT'S a 1D DD
    !IF WE WANT TO SET UP A 2D DD I WILL NEED DO REPLACE Ny_g BY Ny_l
    !AND MAKE OTHER MODIFICATION
    IF(rang == 0)THEN
       S2 = S2 + (overlap-1)
       itN = itN + (overlap-1) * Ny_g
    ELSE IF(rang == Np-1)THEN
       S1 = S1 - (overlap-1)
       it1 = it1 - (overlap-1) * Ny_g
    ELSE
       S1 = S1 - (overlap-1)
       S2 = S2 + (overlap-1)
       it1 = it1 - (overlap-1) * Ny_g
       itN = itN + (overlap-1) * Ny_g
    END IF
   
  End Subroutine CHARGE_OVERLAP
  
end module mod_charge
