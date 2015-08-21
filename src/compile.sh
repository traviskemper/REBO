
rm *~ *.mod
EXE=chosf.x
#ifort main.f -o $EXE 
gfortran  main.f -o $EXE 
##ifort -check all -traceback main.f -o ch.x
