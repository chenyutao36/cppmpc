include ../make.mk

IFLAGS = -I$(EIGEN_HEADER) \
		 -I$(CPPMPC_HEADER) \
	     -I$(QPOASES_HEADER) \
		 -I$(MODEL_HEADER)

MKDIR_P = mkdir -p

OBJ += mpc_common.o
OBJ += casadi_wrapper.o
OBJ += qp_problem.o
OBJ += full_condensing.o
OBJ += qpsolver.o
OBJ += rti_step.o
OBJ += Timer.o

all: directories $(CPPMPC_LIB)/libcppmpc.a $(CPPMPC_LIB)/libcppmpc.so

directories: $(CPPMPC_LIB)

$(CPPMPC_LIB):
	$(MKDIR_P) $(CPPMPC_LIB)
	
$(CPPMPC_LIB)/libcppmpc.so: $(OBJ)
	@echo "Create" $@
	@$(CPP) -shared -o $(CPPMPC_LIB)/libcppmpc.so $(OBJ) $(QPOASES_LINK) $(LINK_DEPENDS)

$(CPPMPC_LIB)/libcppmpc.a: $(OBJ)
	@echo "Create" $@
	@$(AR) $(ARFLAGS) $(CPPMPC_LIB)/libcppmpc.a $(OBJ)

%.o: %.cpp
	@echo "Create" $@
	@$(CPP) $(CPPFLAGS) -c $(IFLAGS) $<

clean: 
	@echo "Cleaning up (src)"
	@rm -rf *.o
