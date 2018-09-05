include ../make.mk

IFLAGS = -I$(MODEL_HEADER)

OBJ = casadi_src.o

all: $(OBJ)

# libmymodel.so: casadi_src.o
# 	@echo "Creating" $@
# 	@$(CPP) -shared -o libmymodel.so casadi_src.o

casadi_src.o:
	@$(CPP) $(CPPFLAGS) -c $(IFLAGS) casadi_src.cpp

clean:
	@echo "Cleaning up (model)"
	@rm -rf *.o *.so