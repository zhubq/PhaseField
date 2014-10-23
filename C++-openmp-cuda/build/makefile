#This is an easier to use and modify makefile, but it is slightly more difficult to read than the simple one:
#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files



#provide the directory where the library Blitz++ installed.
#***********************************
BZDIR=/home/benqiang/blitz-0.10
#***********************************

# define the C++ compiler to use
CXX = g++
# define any compile-time flags
CXXFLAGS =-fopenmp -openmp-report -DBZ_THREADSAFE -ansi -O3 -ftemplate-depth-30

#define the source directory
VPATH=../src
#define the INCLUDE flag.
INCLUDES = -I$(BZDIR)
# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -L$(BZDIR)/lib -L $(VPATH)

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lblitz -lm

#**********TEST IF nvcc IS AVAILABLE
        NVCC_RESULT := $(shell which nvcc 2> NULL)
        
        NVCC_TEST := $(notdir $(NVCC_RESULT))
ifeq ($(NVCC_TEST),nvcc)
        NVCC=nvcc
        LIBS+= -lcudart
        GPUFLAGS= -arch=sm_20
else
        NVCC= g++
        GPUFLAGS=
endif

# define the source files for host compiler
SRCS =main.cpp \
      thermo_data.cpp \
      thermo_functor.cpp \
      phase_field.cpp \
      MicrostructureConstruction.cpp \
      write_data.cpp \
      solver.cpp \
      nucleation.cpp \
      gpu_solver.cpp
# define the object files
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
MAIN = PFM

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
# 
.PHONY: depend clean
# 
all:    $(MAIN)
	@echo   $(MAIN) has been compiled
	@rm -f NULL
	
$(MAIN):$(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJS) -o $(MAIN)  $(LFLAGS) $(LIBS)
# 
# # this is a suffix replacement rule for building .o's from .c's
# # it uses automatic variables $<: the name of the prerequisite of
# # the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# # (see the gnu make manual section about automatic variables)

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS)  $(INCLUDES) $< -o $@
gpu_solver.o:$(VPATH)/gpu_solver.cpp
	@rm -f gpu_solver.o
	$(NVCC) -c $(GPUFLAGS) -O3  $(VPATH)/gpu_solver.cpp
	@echo $(NVCC) is used to compile gpu_solver.cpp
clean:
	$(RM) *.o *~ 
# 
# depend: $(SRCS)
# 	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it