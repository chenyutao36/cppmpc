CC=gcc
Cpp=g++
AR=ar
#CFLAGS=-Wall -O3 -march=native -mfpmath=sse
#CFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__  -D__NO_COPYRIGHT__
CFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -finline-functions -fPIC -DLINUX -D__NO_COPYRIGHT__
ARFLAGS = rcs

QPOASES_PATH = /home/chen/Documents/Packages/qpOASES-3.2.0/include
QPOASES_LIB  = /home/chen/Documents/Packages/qpOASES-3.2.0/bin

QPOASES_LINK = -L$(QPOASES_LIB) -Wl,-rpath=$(QPOASES_LIB) -lqpOASES
LINK_DEPENDS = -llapack -lblas -lm 
CPPMPC_LINK = -L$(CPPMPC_LIB) -Wl,-rpath=$(CPPMPC_LIB) -lcppmpc

CPPMPC_LIB  = /home/chen/Documents/Packages/cppmpc

OBJ += mpc_common.o
OBJ += casadi_src.o
OBJ += casadi_wrapper.o
OBJ += qp_generation.o
OBJ += full_condensing.o
OBJ += qpsolver.o
OBJ += rti_step.o

#OBJ_QPOASES += qpoases_obj/Bounds.o
#OBJ_QPOASES += qpoases_obj/Constraints.o
#OBJ_QPOASES += qpoases_obj/Flipper.o
#OBJ_QPOASES += qpoases_obj/Indexlist.o
#OBJ_QPOASES += qpoases_obj/Matrices.o
#OBJ_QPOASES += qpoases_obj/MessageHandling.o
#OBJ_QPOASES += qpoases_obj/Options.o
#OBJ_QPOASES += qpoases_obj/OQPinterface.o
#OBJ_QPOASES += qpoases_obj/QProblem.o
#OBJ_QPOASES += qpoases_obj/QProblemB.o
#OBJ_QPOASES += qpoases_obj/SolutionAnalysis.o
#OBJ_QPOASES += qpoases_obj/SparseSolver.o
#OBJ_QPOASES += qpoases_obj/SQProblem.o
#OBJ_QPOASES += qpoases_obj/SQProblemSchur.o
#OBJ_QPOASES += qpoases_obj/SubjectTo.o
#OBJ_QPOASES += qpoases_obj/Utils.o

all: lib

# library
#lib: $(OBJ)
#	$(AR) $(ARFLAGS) libcppmpc.a $(OBJ) $(OBJ_QPOASES)

lib: $(OBJ)
	$(Cpp) -shared -o libcppmpc.so $(OBJ) $(QPOASES_LINK) $(LINK_DEPENDS)

qp_generation.o: 
	$(Cpp) $(CFLAGS) -c qp_generation.cpp

casadi_src.o: 
	$(CC) $(CFLAGS) -c casadi_src.c

casadi_wrapper.o: 
	$(Cpp) $(CFLAGS) -c casadi_wrapper.cpp

mpc_common.o: 
	$(Cpp) $(CFLAGS) -c mpc_common.cpp

full_condensing.o: 
	$(Cpp) $(CFLAGS) -c full_condensing.cpp

qpsolver.o: 
	$(Cpp) $(CFLAGS) -c -I$(QPOASES_PATH) -I. qpsolver.cpp

rti_step.o:
	$(Cpp) $(CFLAGS) -c -I$(QPOASES_PATH) -I. rti_step.cpp

#install_static:
#	cp libcppmpc.a /usr/local/lib/libcppmpc.a

#install_shared:
#	cp libcppmpc.so /usr/local/lib/libcppmpc.so
#	cp $(QPOASES_LIB)/libqpOASES.so /usr/local/lib/libqpOASES.so

main: main.o
	$(Cpp) -o main main.o $(CPPMPC_LINK)

main.o:
	$(Cpp) $(CFLAGS) -c -I$(QPOASES_PATH) -I. main.cpp

example1a: example1a.o
	$(Cpp) -o example1a example1a.o $(QPOASES_LINK) $(LINK_DEPENDS)

example1a.o:
	$(Cpp) $(CFLAGS) -c -I$(QPOASES_PATH) example1a.cpp 

clean:
	rm -rf *.o main 

