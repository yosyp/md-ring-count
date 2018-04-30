# The makefile compiles (type make) and makes an executable called mse627
EXEPATH = ..
EXE = $(EXEPATH)/mse627

F90 = ifort
F90FLAGS = -heap-arrays -Ofast -assume buffered_io -ipo
#F90FLAGS = -g -check all -fpe0 -warn -traceback -debug extended
LDFLAGS = 

OBJS = GlobalParameters.o Table.o GlobalVars.o Eftab-pair.o Eftab-SW.o Eftab_TF.o Load_EAM.o Forces.o F-pair.o F-SW.o F_TF.o F_EAM.o Heating.o Nord5.o OpenFiles.o Gather.o Quench.o ReadFiles.o SetInit.o NbList.o Write.o WriteEnd.o WriteInit.o Temper.o Pressure.o Bere-PT.o Vel.o Main.o

# compile and load
default:
	@echo " "
	@echo "Compiling Code MSE627-MD"
	@echo "Version 1.6"
	@echo "FORTRAN 90"
	$(MAKE) $(EXE)

$(EXE):	$(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $(EXE)  $(OBJS)

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90

clean:
	-rm *.o
	-rm *.mod

