include ../make.mk

IFLAGS = -I$(CPPMPC_HEADER) \
	     -I$(QPOASES_HEADER)

OBJ += InvertedPendulum.o

OUT_STATIC += InvertedPendulum_s.exe

OUT_DYNAMIC += InvertedPendulum_d.exe

all: $(OUT_STATIC) $(OUT_DYNAMIC)

$(OUT_DYNAMIC): $(OBJ)
	@echo "Creating" $@
	@$(CPP) -o $@ ../model/casadi_src.o $< $(CPPMPC_LINK)

$(OUT_STATIC): $(OBJ)
	@echo "Creating" $@
	@$(CPP) -o $@ ../model/casadi_src.o $< $(CPPMPC_LIB)/libcppmpc.a $(QPOASES_LIB)/libqpOASES.a -lblas -llapack
	# @$(CPP) -o $@ ../model/casadi_src.o $< -Wl,--start-group $(QPOASES_LIB)/libqpOASES.a $(CPPMPC_LIB)/libcppmpc.a -Wl,--end-group -lblas -llapack

%.o: %.cpp
	@echo "Creating" $@
	@$(CPP) $(CPPFLAGS) -c $(IFLAGS) $<

clean: 
	@echo "Cleaning up (examples)"
	@rm -rf *.o *.exe *.txt
