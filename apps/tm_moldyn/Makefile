#***********************************************************/
# File :  Makefile  (for  moldyn.c)
# Created  :  Shamik Sharma,
# Modified :  Mike Dai Wang 
#***********************************************************
CC = /opt/gcc-4.8/bin/gcc
CXX = /usr/bin/g++

ifeq ($(DEBUG), )
  DEBUG = -O3
endif

#Application options
VERSION = -DSMALL -DTM -DOOP
#    VERSION = -DVERSION_DYNAMIC
#    export VERSION


#   CFLAGS = -O  -dalign -DSMALL -DMEASURE
#   CFLAGS = -O  -dalign -DMEDIUM -DMEASURE
#   CFLAGS = -O  -dalign -DMEDIUM -DPRINT_COORDINATES
#   CFLAGS = -O  -dalign -DMEDIUM -DPRINT_INTERACTION_LIST
#   CFLAGS = -O  -dalign -DMEDIUM 
#   CFLAGS = -O3 -DVERSION_DYNAMIC -DLIBTM -DLARGE -DDEBUG
CXXFLAGS := $(DEBUG) $(VERSION) -mrtm
LDFLAGS = -lm -lpthread -l_tm -lrt 

LIBTM_DIR = ../..
OUT_DIR = .

.c.o:
	$(CXX) $(CXXFLAGS) -c $*.cc 

#***********************************************************
    TARGET = moldyn
    OBJS   = moldyn.o
#***********************************************************

$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(OUT_DIR)/$@ $(OBJS) $(LDFLAGS) -L$(LIBTM_DIR)

#***********************************************************
#    Dependencies
#***********************************************************

moldyn.o   : moldyn.cc Makefile

#***********************************************************
clean:
	/bin/rm -f $(OBJS) $(OUT_DIR)/$(TARGET)
