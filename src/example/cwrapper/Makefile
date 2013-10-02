CC  = gcc
CXX = g++
CFLAGS  = -Wall 
GRACKLE_DEFINES = -DLARGE_INTS  -DCONFIG_BFLOAT_8 
GRACKLE_INCLUDE = -I/home/sleitner/local/include/grackle 
GRACKLE_LIB = -L/home/sleitner/local/lib -lgrackle 

all: example.exe grackle_wrapper.o

grackle_wrapper.o: grackle_wrapper.cxx grackle_wrapper.h
	$(CXX) $(CFLAGS) $(INCLUDES) -I. $(GRACKLE_INCLUDE) -o $@ -c ./$(*F).cxx $(GRACKLE_DEFINES)
example.exe:  grackle_wrapper.o example.c grackle_wrapper.h
	$(CC) $(CFLAGS) -I. $(INCLUDES) grackle_wrapper.o example.c  -o example.exe  $(GRACKLE_LIB) 

clean:
	-@rm -f *.o *.mod *.f *.f90 *~ example.exe
