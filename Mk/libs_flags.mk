OPSYS=$(shell uname)

PLATFORM=$(shell uname -p)
ARCH=.$(PLATFORM)

CXX:=ccache $(CXX)

CXXFLAGS+=-std=c++17

CPPFLAGS+=$(LOCAL_CFLAGS)
LDFLAGS+=$(LOCAL_LDFLAGS)

AGILICIOUS_CFFLAGS=-I./agilicious/agilib/externals/eigen/eigen3 -I./agilicious/agilib/include 
AGILICIOUS_LDFLAGS=-L./agilicious/agilib/build -lagilib -pthread

# CPPFLAGS+=-I./include 
CPPFLAGS+=-I./include $(AGILICIOUS_CFFLAGS) 
LDFLAGS+=$(AGILICIOUS_LDFLAGS) 

#CXXFLAGS+= -g
CXXFLAGS+= -O3 -march=native
#-g -O3 
#-pg
#-march=native 
#CXXFLAGS+= 
#CXXFLAGS+=-std=c++11

