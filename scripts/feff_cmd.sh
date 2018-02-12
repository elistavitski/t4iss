#!/bin/bash

if [ -z "$2" ]
then  wd=`pwd`
    else  wd=$2
fi


if [ -z "$1" ]
then  ncore=1
    else  ncore=$1
fi





case `hostname` in 
    ichost* | icsubmit* )
        module load intel/PSXE2017.u4
        export OMP_NUM_THREADS=1
        FeffPath=~/software/JFEFF/hsw/feff90/linux
#         mpi_cmd="srun -n"
        mpi_cmd="mpirun -np "        
        c=" "
        ;;
    knlhost* | knlsubmit* )
        FeffPath=~/software/JFEFF/knl/feff90/linux
        module load intel/PSXE2017.u4
        export OMP_NUM_THREADS=1        
        mpi_cmd="srun -n"        
        c="-c4"
        ;;
    * )
        FeffPath=~/software/JFEFF/feff90/linux
        mpi_cmd="mpirun -n"        
        c=" "
esac



cd $wd


if [ $ncore -lt 2 ]
then
    > feff.out
    echo "### feff starts at `date`"  >> feff.out     
    echo "### Serial version"         >> feff.out    
    echo "### FeffPath is $FeffPath"  >> feff.out
    echo " "                          >> feff.out    
    $FeffPath/rdinp      >> feff.out
    $FeffPath/atomic     >> feff.out
    $FeffPath/dmdw       >> feff.out
    $FeffPath/pot        >> feff.out
    $FeffPath/opconsat   >> feff.out
    $FeffPath/screen     >> feff.out
    $FeffPath/xsph       >> feff.out
    $FeffPath/fms        >> feff.out
    $FeffPath/mkgtr      >> feff.out
    $FeffPath/path       >> feff.out
    $FeffPath/genfmt     >> feff.out
    $FeffPath/ff2x       >> feff.out
    $FeffPath/sfconv     >> feff.out
    $FeffPath/compton    >> feff.out
    $FeffPath/eels       >> feff.out
#     $FeffPath/ldos       >> feff.out
    echo " "                          >> feff.out 
    echo "### feff ends at `date`"    >> feff.out 
    echo " "                          >> feff.out     



else
    > feff.out
    echo "### feff starts at `date`"  >> feff.out     
    echo "### Parallel version"                   >> feff.out  
    echo "### FeffPath is $FeffPath"              >> feff.out
    echo "### $mpi_cmd $ncore $c \$FeffPath/*"    >> feff.out      
    echo " "                                      >> feff.out    
    $mpi_cmd $ncore $c $FeffPath/MPI/rdinp        >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/atomic       >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/dmdw         >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/pot          >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/opconsat     >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/screen       >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/xsph         >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/fms          >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/mkgtr        >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/path         >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/genfmt       >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/ff2x         >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/sfconv       >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/compton      >> feff.out
    $mpi_cmd $ncore $c $FeffPath/MPI/eels         >> feff.out
#     $mpi_cmd $ncore $c $FeffPath/MPI/ldos         >> feff.out
    echo " "                          >> feff.out 
    echo "### feff ends at `date`"    >> feff.out 
    echo " "                          >> feff.out     
fi






rm -f apot.bin ATOMS.dat band.inp compton.inp config.dat convergence.scf convergence.scf.fine crpa.inp dmdw.inp edges.dat    feff_other_files
rm -f eels.inp emesh.bin emesh.dat feff.bin ff2x.inp fms.bin fms.inp fort.11 fpf0.dat fullspectrum.inp genfmt.inp geom.dat   feff_other_files
rm -f gg.bin gg.dat global.inp gtr.dat hubbard.inp ldos.inp list.dat log1.dat log2.dat log3.dat log4.dat log5.dat log6.dat   feff_other_files
rm -f log.dat logdos.dat logeels.dat logscreen.dat logsfconv.dat mpse.dat opcons.inp paths.dat paths.inp phase.bin pot.bin   feff_other_files
rm -f pot.inp reciprocal.inp rixs.inp screen.inp sfconv.inp sPOSCAR_shifted vtot.dat wscrn.dat xsect.dat xsph.inp .dimensions.dat .feff.error  

rm -f fort.* gtr* ldos* rhoc* atoms.dat


# python -c "import mFEFF; positions, coulomb_matrix, gr, gr_s, xmu = mFEFF.readfeff(plotkey=0)"
