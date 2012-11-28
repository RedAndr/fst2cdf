
RCOMPIL = r.compile
RBUILD  = r.build
FLAGS  = -O 0 

OBJ    = gridder.o gridMercator.o gridPolar.o fst2cdf.o 

LIBS   = rmn_012 netcdf
#LIBPATH  = /local/programs/armnlib/lib

%.o: %.f90
	$(RCOMPIL) -src $< $(FLAGS)

fst2cdf: $(OBJ)
	$(RBUILD) -obj $^ -arch $(ARCH) -abi $(ABI) -opt "=$(LFLAGS)" -o fst2cdf -libappl "$(LIBS)"

clean:
	rm -f fst2cdf *.o *.mod

