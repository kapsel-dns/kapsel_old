#
# $Id: Makefile,v 2.00 2009/08/04 $
#

## default options
ARCH   = linux 
CC     = gcc
CXX    = g++
##CCOPT  = -O2 -static
CCOPT  = -O2 
LINKS  = -lm -lplatform -lstdc++
MKL_PATH = /opt/intel/mkl/10.1.1.019/lib/em64t
MKL_INCLUDE_PATH = /opt/intel/mkl/10.1.1.019/include
GOURMET_HOME_PATH = /usr/local/OCTA2005/GOURMET_2005
GOURMET_LIB_PATH = $(GOURMET_HOME_PATH)/lib/$(ARCH)
GOURMET_INCLUDE_PATH = $(GOURMET_HOME_PATH)/include
OSTYPE = $(shell uname)
#ENV = GCC
ENV = ICC_MKL_OMP
#ENV = MINGW

## options for GCC/CYGWIN
#ifeq ($(ENV), CYGWIN)
ifneq (,$(findstring CYGWIN,$(OSTYPE)))
      ARCH   = cygwin
      CC     = gcc 
      CXX    = g++ 
      CCOPT  = -O3 -fno-inline
      LINKS  = -lm -lplatform 
endif

## options for GCC+MINGW/CYGWIN
ifeq ($(ENV), MINGW)
      ARCH   = win32
      CC     = gcc 
      CXX    = g++ 
      CCOPT  = -O3 -fno-inline -mno-cygwin
      LINKS  = -lm -lplatform 
endif

## options for GCC/LINUX
ifeq ($(ENV), GCC)
      ARCH   = linux 
      CC     = gcc
      CXX    = g++
      CCOPT  = -O3 
      LINKS  = -lm -lplatform -lstdc++
endif

## options for ICC/LINUX
ifeq ($(ENV), ICC)
      ARCH   = linux 
      CC     = icc 
      CXX    = icpc 
      CCOPT  = -O3 -xSSSE3 -axSSSE3 -w0
      LINKS  = -lm -lplatform -lcxaguard -lstdc++
endif

## options for GCC+MKL+OMP/LINUX
ifeq ($(ENV), ICC_MKL_OMP)
      ARCH   = linux 
      CC     = icc 
      CXX    = icpc 
      CCOPT  = -O3 -xSSSE3 -axSSSE3 -ip -openmp -parallel -w0 #-L$(MKL_PATH) -I$(MKL_INCLUDE_PATH) 
      LINKS  = -lplatform -lcxaguard -lstdc++ -lmkl_intel_lp64 -lmkl_intel_thread  -lmkl_core -lm 
endif

CFLAGS 	= $(CCOPT) -L$(GOURMET_LIB_PATH) -I$(GOURMET_INCLUDE_PATH) # -lrfftw -lfftw

OBJS  	= mt19937ar.o\
	operate_electrolyte.o\
	fluct.o\
	alloc.o\
	solute_rhs.o\
	operate_qij.o\
	operate_Qian_Sheng.o\
	fftsg.o\
	fftsg3d.o\
	avs_output.o\
	avs_output_p.o\
	resume.o\
	make_phi.o\
	fluid_solver.o\
	particle_solver.o\
	md_force.o\
	profile.o\
	interaction.o\
	operate_omega.o\
	fft_wrapper.o\
	f_particle.o\
	init_fluid.o\
	init_particle.o\
	input.o\
	sp_3d_ns.o

TARGET 	= kapsel

## Implicit rules

.SUFFIXES: .c .cxx .o .out

## Build rules

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(CFLAGS) $(LINKS)


## Compile

.cxx.o: 
	$(CXX) -c $< $(CFLAGS) -o $@

.c.o: 
	$(CC) -c $< $(CFLAGS) -o $@

## Clean

clean:
	rm -f $(OBJS) $(TARGET) $(TARGET).x
	rm -f *~ *.bak *.x

depend:
	makedepend -- $(CFLAGS) -- *.cxx *.c *.h

