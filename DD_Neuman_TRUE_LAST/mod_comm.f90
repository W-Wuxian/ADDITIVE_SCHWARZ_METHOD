module mod_comm
  use mpi
  use mod_parametres
  Implicit None


CONTAINS

  Subroutine COMM_PTOP(U,UG,UD)
    Real(PR), Dimension(it1:itN), Intent(IN):: U
    Real(PR), Dimension(LBUG:UBUG), Intent(INOUT):: UG
    Real(PR), Dimension(LBUD:UBUD), Intent(INOUT):: UD
!!$    Integer:: A, B, C, D, iD1, iD2, iG1, iG2
!!$    B = itN - 2*(overlap-1)*Ny_g; A = B - (Ny_g-1)
!!$    C = it1 +2*(overlap-1)*Ny_g; D = C + (Ny_g-1)
!!$    iD1 = itN + 1; iD2 = itN + Ny_g
!!$    iG1 = it1 - Ny_g; iG2 = it1 -1
    IF(rang == 0)THEN !                               rang   tag channel
       CALL MPI_SEND(U(A:B),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_RECV(UD(LBUD:UBUD),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
    ELSE IF(rang > 0 .AND. rang<Np-1)THEN
       CALL MPI_RECV(UG(LBUG:UBUG),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
       CALL MPI_SEND(U(C:DD),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_SEND(U(A:B),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_RECV(UD(LBUD:UBUD),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
    ELSE IF(rang == Np-1)THEN
       CALL MPI_RECV(UG(LBUG:UBUG),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
       CALL MPI_SEND(U(C:DD),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
    END IF
  End Subroutine COMM_PTOP


   Subroutine COMM_PTOP_2WAY(U,UG,UD)
    Real(PR), Dimension(it1:itN), Intent(IN):: U
    Real(PR), Dimension(LBUG:UBUG), Intent(INOUT):: UG
    Real(PR), Dimension(LBUD:UBUD), Intent(INOUT):: UD
!!$    Integer:: A, B, C, D, iD1, iD2, iG1, iG2
!!$    B = itN - 2*(overlap-1)*Ny_g; A = B - (Ny_g-1)
!!$    C = it1 +2*(overlap-1)*Ny_g; D = C + (Ny_g-1)
!!$    iD1 = itN + 1; iD2 = itN + Ny_g
!!$    iG1 = it1 - Ny_g; iG2 = it1 -1
     IF(rang == 0)THEN !                               rang   tag channel
       CALL MPI_SEND(U(A:B),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_RECV(UD(LBUD:UBUD),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
    ELSE IF(rang > 0 .AND. rang<Pivot)THEN
       CALL MPI_RECV(UG(LBUG:UBUG),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
       CALL MPI_SEND(U(C:DD),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_SEND(U(A:B),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_RECV(UD(LBUD:UBUD),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
    ELSE IF(rang >= Pivot .AND. rang < Np-1)THEN
       CALL MPI_RECV(UD(LBUD:UBUD),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang,MPI_COMM_WORLD,status,stateinfo)
       CALL MPI_SEND(U(A:B),2*Ny_g,MPI_DOUBLE_PRECISION,rang+1,rang+1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_SEND(U(C:DD),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_RECV(UG(LBUG:UBUG),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
    ELSE IF(rang == Np-1)THEN
       CALL MPI_SEND(U(C:DD),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang-1,MPI_COMM_WORLD,stateinfo)
       CALL MPI_RECV(UG(LBUG:UBUG),2*Ny_g,MPI_DOUBLE_PRECISION,rang-1,rang,MPI_COMM_WORLD,status,stateinfo)
    END IF
  End Subroutine COMM_PTOP_2WAY
  
end module mod_comm
