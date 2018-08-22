CC=gcc
CPP=g++
AR=ar
#CPPFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__  -D__NO_COPYRIGHT__
CPPFLAGS = -Wall -march=native -mfpmath=sse -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -finline-functions -fPIC -DLINUX -D__NO_COPYRIGHT__
ARFLAGS = rcs

QPOASES_PATH = /home/chen/Documents/Packages/qpOASES-3.2.0
CPPMPC_LIB  = /home/chen/Documents/Packages/cppmpc

QPOASES_HEADER = $(QPOASES_PATH)/include
QPOASES_LIB  = $(QPOASES_PATH)/bin
QPOASES_LINK = -L$(QPOASES_LIB) -Wl,-rpath=$(QPOASES_LIB) -lqpOASES
LINK_DEPENDS = -llapack -lblas -lm 
CPPMPC_LINK = -L$(CPPMPC_LIB) -Wl,-rpath=$(CPPMPC_LIB) -lcppmpc

OBJ += mpc_common.o
OBJ += casadi_src.o
OBJ += casadi_wrapper.o
OBJ += qp_generation.o
OBJ += full_condensing.o
OBJ += qpsolver.o
OBJ += rti_step.o
OBJ += Timer.o

all: lib

# library
lib: $(OBJ)
	$(CPP) -shared -o libcppmpc.so $(OBJ) $(QPOASES_LINK) $(LINK_DEPENDS)

qp_generation.o: 
	$(CPP) $(CPPFLAGS) -c qp_generation.cpp

casadi_src.o: 
	$(CC) $(CPPFLAGS) -c casadi_src.c

casadi_wrapper.o: 
	$(CPP) $(CPPFLAGS) -c casadi_wrapper.cpp

mpc_common.o: 
	$(CPP) $(CPPFLAGS) -c mpc_common.cpp

full_condensing.o: 
	$(CPP) $(CPPFLAGS) -c full_condensing.cpp

qpsolver.o: 
	$(CPP) $(CPPFLAGS) -c -I$(QPOASES_HEADER) -I. qpsolver.cpp

rti_step.o:
	$(CPP) $(CPPFLAGS) -c -I$(QPOASES_HEADER) -I. rti_step.cpp

Timer.o:
	$(CPP) $(CPPFLAGS) -c -I. Timer.cpp

main: main.o
	$(CPP) -o main main.o $(CPPMPC_LINK)

main.o:
	$(CPP) $(CPPFLAGS) -c -I$(QPOASES_HEADER) -I. main.cpp

clean:
	rm -rf *.o main

deep_clean:
	rm -rf *.o main *.a *.so

