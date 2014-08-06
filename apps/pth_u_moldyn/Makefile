#***********************************************************/
# File :  Makefile  (for  moldyn.c)
# Created  :  Shamik Sharma,
# Modified :  Mike Dai Wang 
#***********************************************************
CC = /opt/gcc-4.8/bin/gcc
CXX = /opt/gcc-4.8/bin/g++

ifeq ($(DEBUG), )
  DEBUG = -O3
endif

#Application options
VERSION = -DSUPER
#    VERSION = -DVERSION_DYNAMIC
#    export VERSION


#   CFLAGS = -O  -dalign -DSMALL -DMEASURE
#   CFLAGS = -O  -dalign -DMEDIUM -DMEASURE
#   CFLAGS = -O  -dalign -DMEDIUM -DPRINT_COORDINATES
#   CFLAGS = -O  -dalign -DMEDIUM -DPRINT_INTERACTION_LIST
#   CFLAGS = -O  -dalign -DMEDIUM 
#   CFLAGS = -O3 -DVERSION_DYNAMIC -DLIBTM -DLARGE -DDEBUG
CXXFLAGS := $(DEBUG) $(VERSION)
LDFLAGS = -lm -lpthread -lrt -l_tm

LIBTM_DIR = ../..
OUT_DIR = ../.. 

.c.o:
	$(CXX) $(CXXFLAGS) -c $*.cc 

#***********************************************************
    TARGET = moldyn
    OBJS   = moldyn.o
#***********************************************************

$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) -L$(LIBTM_DIR)

#***********************************************************
#    Dependencies
#***********************************************************

moldyn.o   : moldyn.cc Makefile

#***********************************************************
clean:
	/bin/rm -f $(OBJS) $(TARGET)
