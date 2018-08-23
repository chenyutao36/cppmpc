include ../make.mk

IFLAGS = -I$(CPPMPC_HEADER) \
	     -I$(QPOASES_HEADER) \
		 -I$(MODEL_HEADER)

OBJ += mpc_common.o
OBJ += casadi_wrapper.o
OBJ += qp_generation.o
OBJ += full_condensing.o
OBJ += qpsolver.o
OBJ += rti_step.o
OBJ += Timer.o

all: $(CPPMPC_LIB)/libcppmpc.so
	
$(CPPMPC_LIB)/libcppmpc.so: $(OBJ)
	@echo "Create" $@
	@$(CPP) -shared -o $(CPPMPC_LIB)/libcppmpc.so $(OBJ) $(QPOASES_LINK) $(LINK_DEPENDS)

%.o: %.cpp
	@echo "Create" $@
	@$(CPP) $(CPPFLAGS) -c $(IFLAGS) $<

clean: 
	@echo "Cleaning up (src)"
	@rm -rf *.o
