CC=gcc
Cpp=g++
AR=ar
CFLAGS=-Wall -fPIC -O#-march=native -mfpmath=sse
ARFLAGS = crs

QPOASES_PATH = /home/chenyutao/Documents/Packages/qpOASES-3.2.1/include
# QPOASES_LIB  = /home/chenyutao/Documents/Packages/qpOASES-3.2.1/bin
QPOASES_LIB = .

OBJ_main = main.o

OBJ += mpc_common.o
OBJ += casadi_src.o
OBJ += casadi_wrapper.o
OBJ += qp_generation.o
OBJ += full_condensing.o
OBJ += qpsolver.o
OBJ += rti_step.o


all: lib main

# for library

# lib: $(OBJ)
# 	$(AR) $(ARFLAGS) libcppmpc.a $(OBJ) 
lib: $(OBJ)
	$(Cpp) -shared $(OBJ) -L$(QPOASES_LIB) -o libcppmpc.so -lqpOASES 

qp_generation.o: qp_generation.cpp qp_generation.hpp
	$(Cpp) $(CFLAGS) -c qp_generation.cpp

casadi_src.o: casadi_src.c casadi_src.h 
	$(CC) $(CFLAGS) -c casadi_src.c

casadi_wrapper.o: casadi_wrapper.cpp casadi_wrapper.hpp
	$(Cpp) $(CFLAGS) -c casadi_wrapper.cpp

mpc_common.o: mpc_common.cpp mpc_common.hpp
	$(Cpp) $(CFLAGS) -c mpc_common.cpp

full_condensing.o: full_condensing.cpp full_condensing.hpp
	$(Cpp) $(CFLAGS) -c full_condensing.cpp

qpsolver.o: qpsolver.cpp qpsolver.hpp
	$(Cpp) $(CFLAGS) -c -I$(QPOASES_PATH) qpsolver.cpp
#-L$(QPOASES_LIB) -lqpOASES

rti_step.o: rti_step.cpp rti_step.hpp
	$(Cpp) $(CFLAGS) -c rti_step.cpp

# for main

main: $(OBJ_main)
	$(Cpp) -L. -o main $(OBJ_main) -lqpOASES -lcppmpc 

main.o: main.cpp
	$(Cpp) -Wall -O -c main.cpp

clean:
	rm -rf *.o main 

