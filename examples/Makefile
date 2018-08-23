include ../make.mk

IFLAGS = -I$(CPPMPC_HEADER) \
	     -I$(QPOASES_HEADER)

OBJ += InvertedPendulum.o

OUT += InvertedPendulum.exe

all: $(OUT)

$(OUT): $(OBJ)
	@echo "Creating" $@
	@$(CPP) -o $@ $< $(CPPMPC_LINK) $(MODEL_LINK)

%.o: %.cpp
	@echo "Creating" $@
	@$(CPP) $(CPPFLAGS) -c $(IFLAGS) $<

clean: 
	@echo "Cleaning up (examples)"
	@rm -rf *.o *.exe
