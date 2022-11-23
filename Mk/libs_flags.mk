OPSYS=$(shell uname)

PLATFORM=$(shell uname -p)
ARCH=.$(PLATFORM)

CXX:=ccache $(CXX)

CXXFLAGS+=-std=c++17

CPPFLAGS+=$(LOCAL_CFLAGS)
LDFLAGS+=$(LOCAL_LDFLAGS)

#AGILICIOUS_CFFLAGS=-I./agilicious/agilib/externals/eigen/eigen3 -I./agilicious/agilib/include
#AGILICIOUS_LDFLAGS=-L./agilicious/agilib/build -lgtest -lpthread -lagilib -lyaml-cpp -leigen

# CPPFLAGS+=-I./include 
CPPFLAGS+=-I./include -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I/opt/homebrew/Cellar/yaml-cpp/0.7.0/include
LDFLAGS+=-lpthread -L/opt/homebrew/Cellar/yaml-cpp/0.7.0/lib -lyaml-cpp

#CXXFLAGS+= -g
CXXFLAGS+= -O3
#-g -O3 
#-pg
#-march=native 
#CXXFLAGS+= 
#CXXFLAGS+=-std=c++11

