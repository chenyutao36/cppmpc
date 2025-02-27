CC=gcc
CPP=g++
AR=ar
CPPFLAGS = -Wall -O3 -march=native -mfpmath=sse -fPIC -std=c++11
ARFLAGS = rcs

QPOASES_PATH = /home/cyt/Documents/qpOASES-3.2.0
EIGEN_PATH = /home/cyt/Documents/eigen-3.4.0
LINALG_PATH = /usr/lib/x86_64-linux-gnu
CPPMPC_PATH  = $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

EIGEN_HEADER = $(EIGEN_PATH)
QPOASES_HEADER = $(QPOASES_PATH)/include
QPOASES_LIB  = $(QPOASES_PATH)/bin
CPPMPC_HEADER = $(CPPMPC_PATH)/include
CPPMPC_LIB  = $(CPPMPC_PATH)/lib
QPOASES_LINK = -L$(QPOASES_LIB) -Wl,-rpath=$(QPOASES_LIB) -lqpOASES
LINK_DEPENDS = -L$(LINALG_PATH)/ -llapack -lblas -lm 
CPPMPC_LINK = -L$(CPPMPC_LIB) -Wl,-rpath=$(CPPMPC_LIB) -lcppmpc
MODEL_HEADER = $(CPPMPC_PATH)/model
MODEL_LINK = -L$(CPPMPC_PATH)/model -Wl,-rpath=$(CPPMPC_PATH)/model -lmymodel

