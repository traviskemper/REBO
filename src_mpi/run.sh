#!/bin/bash
#
#$ -N tPMMAv4
#$ -cwd
#$ -S /bin/bash
#$ -e error.out
#$ -o screen.out
#$ -pe mpi 2
#$ -q debug.q
# -q production.q
# -q highmem.q

##Set variables for run
HD=/home/travisk/CHOSp1/PMMAv4
SDIR=/state/partition1/$JOB_NAME    #Scratch directory on node
RID=t1
EXE=main.x                       #Executable REBO program
LOGFILE=${JOB_NAME}.log
export RUNPATH=$HD

## Create subdirectory
mkdir $HD/$RID

## Record the 
##    node ID, Job name, Job number and date in log file
echo "Master node:" $HOSTNAME >> $RID/$LOGFILE
echo  "Job name:" $JOB_NAME >> $RID/$LOGFILE
echo  "Job number:" $JOB_ID  >> $RID/$LOGFILE
echo "Start time:" `date`  >> $RID/$LOGFILE
echo ""


#Create Scratch directory
##mkdir -p $SDIR/Execute

##Copy over need files for a REBO run
##cp -r ../Commons $SDIR/Commons
##cp -r ../Spline $SDIR/Spline
##cp $HD/$EXE $SDIR/Execute/
##cp $HD/* $SDIR/Execute/
#
##Change to scratch directory
##cd $SDIR/Execute

#Exicute REBO
/opt/intel/mpich-1.2.7p1/bin/mpirun -machinefile $TMPDIR/machines -np $NSLOTS $RUNPATH/$EXE &> $JOB_NAME.out


#Copy the scratch files back to the home directory and remove the scratch directory from the compute node:
##cp -ru  $SDIR/Execute/* $HD/$RID/
##cd $HD
##rm fort.*  !Remove any unwanted files from exe

## Record finish time
echo "End time:" `date`  >> $RID/$LOGFILE

## !!Remove strach files otherwise node will fill up and be become useless!!
##rm -rf $SDIR

