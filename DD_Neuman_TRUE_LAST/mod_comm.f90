module mod_comm
  use mpi
  use mod_parametres
  Implicit None


CONTAINS

  Subroutine COMM_PTOP(U,UG,UD)
    Real(PR), Dimension(it1:itN), Intent(IN):: U
    Real(PR), Dimension(LBUG:UBUG), Intent(INOUT):: UG
    Real(PR), Dimension(LBUD:UBUD), Intent(INOUT):: UD

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
