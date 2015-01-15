
# Makefile for clang/libc++ projects

#******************************************************************************

LINK=llvm-link
CC=clang
CXX=clang++
AR=llvm-ar
AS=llvm-as
NM=llvm-nm

#******************************************************************************

INCLUDES += -I../src
INCLUDES += -I../src/yson/src

VPATH = ./:../src

DEFINES = 

# Avoid errors in math.h: "unknown type name '__extern_always_inline'"
DEFINES += -D__extern_always_inline=inline

# INCLUDES += -I/llvm-svn/include/c++/v1
# LIBS += -L/llvm-svn/lib
LIBS += -lc++


# CFLAGS = -stdlib=libc++ -g -O3 -ffast-math -msse3
CFLAGS = -stdlib=libc++ -g -O3 -ffast-math
CFLAGS += $(DEFINES) $(INCLUDES)
CFLAGS += -Wno-deprecated-register

PWD = $(shell pwd)

ifeq ($(shell uname -s), Darwin)
	CFLAGS += -DMACOSX
	LIBS += -framework OpenGL -framework Cocoa
else
	CFLAGS += -DLINUX
	LIBS += -lGL -lGLU
endif

CXXFLAGS = -std=c++11 $(CFLAGS)


#******************************************************************************
# Dependency rules
#******************************************************************************

.PHONY: all default clean depend echo none disasm

default: test

install:

clean:
	rm -f basic_tests
	rm -f allocator_tests
	rm -f quadedge_tests

test: allocator_tests basic_tests

basic_tests: basic_tests.cpp lambdatest.h
	$(CXX) $(CXXFLAGS) basic_tests.cpp -o basic_tests
	./basic_tests

quadedge_tests: quadedge_tests.cpp quadedge.h lambdatest.h
	$(CXX) $(CXXFLAGS) quadedge_tests.cpp -o quadedge_tests
	./quadedge_tests

trimesh_tests: trimesh_tests.cpp ../src/trimesh.cpp trimesh.h lambdatest.h
	$(CXX) $(CXXFLAGS) trimesh_tests.cpp ../src/trimesh.cpp -o trimesh_tests
	./trimesh_tests


utilities_tests: ../src/utilities.h utilities_tests.cpp lambdatest.h
	$(CXX) $(CXXFLAGS) utilities_tests.cpp -o utilities_tests
	./utilities_tests


#******************************************************************************
# End of file
#******************************************************************************
