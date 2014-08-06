####################################################################################################
##
## CONFIGURATION
##
####################################################################################################

## Compilation parameters
CC = /opt/gcc-4.8/bin/gcc
CXX = /opt/gcc-4.8/bin/g++

ifeq ($(DEBUG), )
	DEBUG = -O3
endif

#Application options
VERSION = -DTM -DTIMER

#LibTM Options
OPTIONS = -D_MAX_RETRIES=0

CXXFLAGS := -Wall $(DEBUG) $(VERSION) $(OPTIONS) -mrtm
LDFLAGS = -lpthread -lstdc++ -l_tm


## Source and binary files
SRCS      =   $(wildcard *.cc)
OBJS      =   $(patsubst %.cc,%.o,$(SRCS))

LIBTM_DIR =		../..
OUT_DIR		=		../..
TARGET		=		array

$(OBJS)		:	$(SRCS) Makefile

$(TARGET)	:	$(OBJS)
						$(CXX) $(CFLAGS) -o $(OUT_DIR)/$@ $(OBJS) $(LDFLAGS) -L$(LIBTM_DIR)

default		: $(TARGET)

clean:
		-rm -f $(OBJS) $(OUT_DIR)/$(TARGET)

.PHONY: clean
