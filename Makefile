FC      = gfortran
FCFLAGS = -O3 -fopenmp
LIBS = -lopenblas -lfftw3

vpath %.f90 .

all: main

main: algor.o cell.o checkpoint.o constants.o eam5.o fourier_interpolation.o fourier_pes.o hash.o geometry.o io.o lj.o md.o pathint.o pes_scan.o phonon.o potential.o singlepoint.o

algor.o: constants.o hash.o io.o

cell.o: constants.o io.o algor.o checkpoint.o

checkpoint.o: io.o

constants.o: io.o

eam5.o: constants.o cell.o algor.o io.o

fourier_interpolation.o: constants.o io.o

fourier_pes.o: constants.o cell.o fourier_interpolation.o io.o

geometry.o: constants.o cell.o io.o potential.o

hash.o:

io.o: hash.o

lj.o: constants.o io.o cell.o

md.o: constants.o io.o cell.o potential.o algor.o checkpoint.o

pathint.o: constants.o cell.o md.o io.o algor.o checkpoint.o

pes_scan.o: constants.o cell.o io.o potential.o

phonon.o: constants.o cell.o io.o potential.o

potential.o: constants.o cell.o lj.o eam5.o fourier_pes.o io.o

singlepoint.o: constants.o io.o cell.o potential.o geometry.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LIBS)

clea: clean

clena: clean

clean:
	rm -rf *.o *.mod main
