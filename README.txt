REBO codes:

The src and src_mpi contain the code for the serial and parallel version respectively. Each can be compiled using
 ./compile.sh
as long as ifort or gfortran is in your PATH.

2) Running 
Both codes need the Spline directory to be located ../ 
The code use an input file, which contains KEYWORD (space) = (space) VALUE and can have any name.
The code can be ran using 
./chosf-s.x < input.d 

3) Bugs


