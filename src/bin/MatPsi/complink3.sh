#!/bin/sh 

# ABOUT PSI4 LIBRARIES  
# 1. All $psi4 mentioned below better be /home/haichen/psi4_for_MatPsi/ 
#    as the original Psi4 source code has been changed a little bit. 
# 2. After compiling Psi4, go to $psi4/objdir/include/libint, recompile 
#    vrr_build.c by adding a -fPIC flag, then replace the original vrr_build.o
#    in $psi4/objdir/lib/libPSI_int.a with the one newly compiled by ar r. 
#    Refer to $psi4/objdir/out about how to compile vrr_build.c. This allows
#    Mex to compile our code into a dynamic library file required by Matlab. 

# COMPILATION NOTICE 
# Please make sure: 
# (1) ./share contains read_options.o from $psi4/objdir/src/bin/psi4/ 
# (2) ./share/lib contains all the libXXX folders from $psi4/src/lib 
#       and all the libPSI_XXX.a static libraries from $psi4/objdir/lib/ 
# (3) ./share/include contains all the .h files from $psi4/include 
#       and folders libderiv and libint from $psi4/objdir/include 
# (4) ./share/basis contains all the .gbs basis set files from $psi4/lib/basis

# Simply execute this file to compile (and link). 
# Refer to original psi4 documents if met blas/lapack/boost issues. 

mex -largeArrayDims -DLinux -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE CFLAGS='$CFLAGS -O2 -fopenmp -fPIC -pthread' MatPsi_mex.cpp -I./share/lib -I/usr/include/python2.7 -I./share/include -I/usr/include -I/home/haichen/program/lib/boost/include CXXFLAGS='$CXXFLAGS -O2 -fopenmp -fPIC -pthread' LDFLAGS='$LDFLAGS -fopenmp -fPIC -pthread -Wl,-export-dynamic -Wl,--whole-archive ./share/read_options.o ./share/lib/libPSI_scf_solver.a ./share/lib/libPSI_fock.a ./share/lib/libPSI_functional.a ./share/lib/libPSI_3index.a ./share/lib/libPSI_mints.a ./share/lib/libPSI_trans.a ./share/lib/libPSI_dpd.a ./share/lib/libPSI_chkpt.a ./share/lib/libPSI_iwl.a ./share/lib/libPSI_psio.a ./share/lib/libPSI_qt.a ./share/lib/libPSI_ciomr.a ./share/lib/libPSI_options.a ./share/lib/libPSI_util.a ./share/lib/libPSI_deriv.a ./share/lib/libPSI_int.a ./share/lib/libPSI_diis.a ./share/lib/libPSI_parallel.a -Wl,--no-whole-archive -lm -llapack -lcblas -lf77blas -latlas -L/home/haichen/program/lib/atlas/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7/../../.. -lgfortran -lm -lf77blas -latlas -L/home/haichen/program/lib/atlas/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.7/../../.. -lgfortran -lm /home/haichen/program/lib/boost/lib/libboost_filesystem.a /home/haichen/program/lib/boost/lib/libboost_python.a /home/haichen/program/lib/boost/lib/libboost_regex.a /home/haichen/program/lib/boost/lib/libboost_serialization.a /home/haichen/program/lib/boost/lib/libboost_system.a /home/haichen/program/lib/boost/lib/libboost_thread.a -lrt -lpthread -L/usr/lib/python2.7/config-x86_64-linux-gnu -L/usr/lib -lpthread -ldl -lutil -lm -lpython2.7 -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions' -cxx 


