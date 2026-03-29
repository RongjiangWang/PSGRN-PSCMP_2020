FC = mpif90
FFLAGS ?= -O3 -fopenmp
TARGET := ../bin/psgrn2020_mpi
SOURCES := $(sort $(wildcard *.f))

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES) | ../bin
	$(FC) $(FFLAGS) $(SOURCES) -o $@

../bin:
	mkdir -p $@

clean:
	rm -f $(TARGET) *.o *.mod
