#########################################################
# REPLACE the main_prg_name below with the name of
# your main program without the ".f90" extension
# REPLACE the module_name below with your module
# name without the ".f90" extension
# PROGRAM NAME (no .f90)
BIN=./bin
OBJ=./obj
VPATH=./src
main=PFM
# MODULE NAME (no .f90)
mod1=sparse
mod2=pre_process
mod3=module_nucleation
mod4=util_func
#mod5=post_process
mod6=common_variable
#########################################################-O3 -ipo -xHost 
cmplr=ifort
src_files=$(wildcard $(VPATH)/*.f90)
dir=$(notdir $(src_files))
obj_files=$(patsubst %.f90,$(OBJ)/%.o,$(dir))

option=-fopenmp -ffree-line-length-200 -fdefault-real-8 -fdefault-double-8 -O3 ### gfortran compilaion options
option=-openmp -fast -xHost -autodouble -module $(OBJ) 
##-g #intel fortran options
$(BIN)/$(main):		$(obj_files)
		$(cmplr) $(option) $(obj_files) -o $(BIN)/$(main)
$(OBJ)/$(main).o:	$(main).f90 $(OBJ)/$(mod1).o $(OBJ)/$(mod2).o $(OBJ)/$(mod3).o $(OBJ)/$(mod4).o $(OBJ)/$(mod6).o
		$(cmplr) $(option) -c $< -o $@ -I $(OBJ)
$(OBJ)/$(mod1).o:	$(mod1).f90 $(OBJ)/$(mod4).o
		$(cmplr) $(option) -c $< -o $@ -I $(OBJ)
$(OBJ)/$(mod2).o:	$(mod2).f90 $(OBJ)/$(mod4).o
		$(cmplr) $(option) -c $< -o $@ -I $(OBJ)	
$(OBJ)/$(mod3).o:	$(mod3).f90 $(OBJ)/$(mod4).o
		$(cmplr) $(option) -c $< -o $(OBJ)/$(mod3).o -I $(OBJ)
$(OBJ)/$(mod4).o:	$(mod4).f90 $(OBJ)/$(mod6).o
		$(cmplr) $(option) -c $< -o $@ -I $(OBJ)	
$(OBJ)/$(mod6).o:	$(mod6).f90
		$(cmplr) $(option) -c $< -o $@


clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod