#
# $Id: Makefile,v 3.00 2013/01/13 $
#

## default options
ARCH   = linux 
AUX= ./Tools
GOURMET_HOME_PATH = /usr/local/OCTA2010/GOURMET_2010
GOURMET_LIB_PATH = $(GOURMET_HOME_PATH)/lib/$(ARCH)
GOURMET_INCLUDE_PATH = $(GOURMET_HOME_PATH)/include

OSTYPE = $(shell uname)
#ENV = GCC
#ENV = ICC
#ENV = ICC_MKL_OMP
#ENV = MINGW
#ENV = CYGWIN

## options for GCC/CYGWIN
ifeq ($(ENV), CYGWIN)
#ifneq (,$(findstring CYGWIN,$(OSTYPE)))
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
      ARCH   = linux_64
      CC     = gcc
      CXX    = g++
      CCOPT  = -O3 
      LINKS  = -lm -lplatform -lstdc++

      GOURMET_LIB_PATH = $(GOURMET_HOME_PATH)/lib/$(ARCH)_$(CC)
      GOURMET_INCLUDE_PATH = $(GOURMET_HOME_PATH)/include
endif

## options for ICC/LINUX
ifeq ($(ENV), ICC)
      ARCH   = linux_64
      CC     = icc 
      CXX    = icpc 
      CCOPT  = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -w0
      LINKS  = -lm -lplatform -lcxaguard -lstdc++

      GOURMET_LIB_PATH = $(GOURMET_HOME_PATH)/lib/$(ARCH)_$(CC)
      GOURMET_INCLUDE_PATH = $(GOURMET_HOME_PATH)/include
endif

## options for GCC+MKL+OMP/LINUX
ifeq ($(ENV), ICC_MKL_OMP)
      MKL_DIR = /home/opt/intel/composer_xe_2013/mkl
      MKL_PATH= $(MKL_DIR)/lib/intel64
      MKL_INCLUDE_PATH = $(MKL_DIR)/include

      ARCH   = linux_64
      CC     = icc 
      CXX    = icpc 
      CCOPT  = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2\
	-ip -openmp -parallel -w0 -L$(MKL_PATH) -I$(MKL_INCLUDE_PATH) 
      LINKS  = -lplatform -lcxaguard -lstdc++\
	-lmkl_intel_lp64 -lmkl_intel_thread  -lmkl_core -lm
      GOURMET_LIB_PATH = $(GOURMET_HOME_PATH)/lib/$(ARCH)_$(CC)
      GOURMET_INCLUDE_PATH = $(GOURMET_HOME_PATH)/include
endif

CFLAGS 	= $(CCOPT) -L$(GOURMET_LIB_PATH) -I$(GOURMET_INCLUDE_PATH)

OBJS  	= mt19937ar.o\
	operate_electrolyte.o\
	fluct.o\
	alloc.o\
	solute_rhs.o\
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
	rigid_body.o\
	operate_surface.o\
	sp_3d_ns.o

XYZ_OBJS= alloc.o\
	rigid_body.o\
	$(AUX)/udf2xyz.o

TARGET 	= kapsel
XYZ	= udf2xyz

## Implicit rules

.SUFFIXES: .c .cxx .o .out

## Build rules

all: $(TARGET) $(XYZ)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(CFLAGS) $(LINKS)

$(XYZ): $(XYZ_OBJS)
	$(CXX) $(XYZ_OBJS) -o $(XYZ) $(CFLAGS) $(LINKS)


## Compile

.cxx.o: 
	$(CXX) -c $< $(CFLAGS) -o $@

.c.o: 
	$(CC) -c $< $(CFLAGS) -o $@

## Clean

clean:
	rm -f $(OBJS) $(AUX)/$(XYZ_OBJS) $(TARGET) $(XYZ)
	rm -f *~ *.bak

depend:
	makedepend -- $(CFLAGS) -- *.cxx *.c *.h
