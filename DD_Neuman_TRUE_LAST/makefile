#F90 = gfortran
F90 = mpif90
FLAG = -O2 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

EXE = run

all : $(EXE)

$(EXE) : mod_parametres.o mod_charge.o mod_comm.o mod_sol.o mod_fonctions.o mod_gradconj.o mod_gmres.o mod_constr_mat.o mod_matsys.o main.o
	$(F90) $(FLAG) -o $(EXE) $^

%.o : %.f90
	$(F90) $(FLAG) -c $<

clean :
	rm -f *.o *.mod fort.* *.dat

cleanREP :
	rm -r ./MAT_VIZU/*.dat ./CHARGE/*.dat ./SOL_NUMERIQUE/*.dat ./SOL_NUMERIQUE/U* ./EXACTE_ERREUR/*.dat *~
