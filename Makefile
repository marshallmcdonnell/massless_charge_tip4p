# compiler
FC = gfortran

# compile flags
FCFLAGS = -g 
# link flags
FLFLAGS = -g 

# source files and objects
SRC = $(patsubst %.f, %.o, $(wildcard *.f)) 

# program name
PROGRAM = add_massless_charges 

all: $(PROGRAM)

$(PROGRAM): $(SRC)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f
	$(FC) $(FLFLAGS) -c $^

clean:
	rm -f *.o *.mod
