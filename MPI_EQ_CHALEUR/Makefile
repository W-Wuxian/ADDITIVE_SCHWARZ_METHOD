
#paramètres

F90=mpif90
##FFLAGS=-o -wall -ffpe-trap=invalid -pedantic -g
#ou
DEBUG=-O3 #-wall -ffpe-trap=invalid -fbounds-check -pedantic  -g

COMPILE= mpif90 $(DEBUG)

#edition de liens creation executable nommé run
run : main.o ModData.o ModF.o ModGC.o ModSOL.o
	$(COMPILE)-o run main.o ModData.o ModF.o ModGC.o ModSOL.o
	@echo compilation terminée

#compilation
main.o : main.f90 ModData.o ModF.o ModGC.o ModSOL.o
	$(COMPILE) -c main.f90
ModData.o : ModData.f90
	$(COMPILE) -c ModData.f90
ModF.o : ModF.f90
	$(COMPILE) -c ModF.f90
ModGC.o : ModGC.f90
	$(COMPILE) -c ModGC.f90
ModSOL.o : ModSOL.f90
	$(COMPILE) -c ModSOL.f90
#destrution des fichiers objets et modules
clean:
	rm -f -- -wall ./-wall.o *.o *.mod *.dat *.txt *.vtk  *~
cleanREP:
	rm -r ./M1/*.dat ./F/*.dat ./EXACTE_ERREUR/*.dat ./SOL_NUMERIQUE/*


#gnufor2.o : gnufor2.f90
#	$(COMPILE) -c gnufor2.f90
#	$(F90) -c gnufor2.f90


