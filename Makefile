
include make.mk


all: src model examples

src:
	@cd $@; ${MAKE} -s

model:
	@cd $@; ${MAKE} -s

examples: src
	@cd $@; ${MAKE} -s

clean:
	@cd src      && ${MAKE} -s clean
	@cd examples && ${MAKE} -s clean
	@cd lib      && echo "Cleaning up (libraries)" && rm -rf *.so
	@echo "Cleaning up (data)" && rm -rf *.txt
	@cd model    && ${MAKE} -s clean

.PHONY: all src model examples
