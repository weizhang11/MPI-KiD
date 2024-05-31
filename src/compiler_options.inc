

VPATH = ./:./UM_source

#----------------------------------------------------------------------#
### cray ###############################################################
#----------------------------------------------------------------------#
ifeq ($(VENDOR),cray)
NCINC=/opt/cray/pe/netcdf/4.8.1.3/CRAYCLANG/14.0/include
NCLIB=/opt/cray/pe/netcdf/4.8.1.3/CRAYCLANG/14.0/lib
MPIINC=/opt/cray/pe/mpich/8.1.16/ofi/cray/10.0/include
MPILIB=/opt/cray/pe/mpich/8.1.16/ofi/cray/10.0/lib

ifeq ($(DEF_DPREC), true)
FPP = -s real64 -eZ -DPI=pi -h pic -h list=a   #CRAY
else
FPP = -eZ -DPI=pi -h pic -h list=a   #CRAY
endif

endif
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
### nvidia #############################################################
#----------------------------------------------------------------------#
ifeq ($(VENDOR),nvidia)
NCINC=/opt/cray/pe/netcdf/4.8.1.3/NVIDIA/20.7/include
NCLIB=/opt/cray/pe/netcdf/4.8.1.3/NVIDIA/20.7/lib
MPIINC=/opt/cray/pe/mpich/8.1.16/ofi/nvidia/20.7/include
MPILIB=/opt/cray/pe/mpich/8.1.16/ofi/nvidia/20.7/lib
TESTING=

ifeq ($(DEF_DPREC), true)
FPP = -r8 -Mpreprocess -DPI=pi -mcmodel=medium -Mlarge_arrays  #nvidia
else
FPP = -Mpreprocess -DPI=pi -mcmodel=medium -Mlarge_arrays  #nvidia
endif

endif
#----------------------------------------------------------------------#
#----------------------------------------------------------------------#
### nvidia gpu ########################################################
#----------------------------------------------------------------------#
ifeq ($(VENDOR),nvidia_gpu)

NCINC=/opt/cray/pe/netcdf/4.9.0.3/NVIDIA/20.7/include
NCLIB=/opt/cray/pe/netcdf/4.9.0.3/NVIDIA/20.7/lib
MPIINC=/opt/cray/pe/mpich/8.1.24/ofi/nvidia/20.7/include
MPILIB=/opt/cray/pe/mpich/8.1.24/ofi/nvidia/20.7/lib

ifeq ($(DEF_DPREC), true)
FPP = -r8 -Mpreprocess -DPI=pi -mcmodel=medium -acc -ta=tesla:managed -Minfo=accel -Mlarge_arrays  
#FPP = -r8 -Mpreprocess -DPI=pi -mcmodel=medium -acc=gpu -Minfo=accel -Mlarge_arrays 
else
FPP = -Mpreprocess -DPI=pi -mcmodel=medium -acc=gpu -Minfo=accel -Mlarge_array
endif

endif

#----------------------------------------------------------------------#

#----------------------------------------------------------------------#
### gnu (need to configure) ############################################
#----------------------------------------------------------------------#
ifeq ($(VENDOR),gnu)
#NCINC=/opt/cray/pe/netcdf/4.8.1.3/NVIDIA/20.7/include
#NCLIB=/opt/cray/pe/netcdf/4.8.1.3/NVIDIA/20.7/lib
#MPIINC=/opt/cray/pe/mpich/8.1.16/ofi/nvidia/20.7/include
#MPILIB=/opt/cray/pe/mpich/8.1.16/ofi/nvidia/20.7/lib
TESTING=
#FPP = -DPI=pi -cpp #GNU
endif
#----------------------------------------------------------------------#

#ifeq ($(COMPILER), ftn)
#endif

OPTS = -O3 -g 
#OPTS = -O0 -g

ifdef PCOMPILER
OPTS+= -DUSE_MPI 
F90=$(PCOMPILER) $(OPTS) 
LD=$(F90) -acc # Linker

INCDF=-I$(NCINC) -I$(MPIINC)
### mpif for SUMMIT
#LNCDF=-L$(NCLIB) -lnetcdff -I$(MPILIB) -lmpi_ibm_mpifh
### mpif for MILLER
LNCDF=-L$(NCLIB) -lnetcdff -I$(MPILIB) -lmpichf90

else

F90=$(COMPILER) $(OPTS)
LD=$(F90)  # Linker
INCDF=-I$(NCINC)
LNCDF=-L$(NCLIB) -lnetcdff
endif

LDFLAGS=$(INCDF) $(LNCDF)

EXECUTABLE=$(MODEL)_$(CASE).exe

F2PY=/usr/local/bin/f2py 
F2PY_OPTS=--fcompiler-exec=$(COMPILER) $(LDFLAGS) 
