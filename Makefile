F90    = gfortran
FFLAGS = -O3 -ffpe-trap=invalid,zero,overflow -ffree-line-length-none

OBJS = opkda1.o opkda2.o opkdmain.o chemistry_mod.o aerosol_mod.o main.o 
EXE  = main.exe

all: $(EXE)

$(EXE): $(OBJS)
	$(F90) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

%.o: %.f
	$(F90) $(FFLAGS) -std=legacy -w -c $< -o $@

clean:
	@rm -vf *.o *.mod *.exe
