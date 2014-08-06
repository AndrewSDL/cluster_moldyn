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
VERSION = -DTM

#LibTM Options
OPTIONS = 

CXXFLAGS := -Wall $(DEBUG) $(VERSION) $(OPTIONS) -mrtm -D_GNU_SOURCE
LDFLAGS = -pthread -lstdc++


## Source and binary files
BIN_DIR   =   build
SRC_DIR   =   src
APP_DIR   =   apps

C_SRCS    =   $(wildcard $(SRC_DIR)/*.c)
C_OBJS    =   $(patsubst $(SRC_DIR)/%.c,$(BIN_DIR)/%.o,$(C_SRCS)) 
CC_SRCS   =   $(wildcard $(SRC_DIR)/*.cc)
CC_OBJS   =   $(patsubst $(SRC_DIR)/%.cc,$(BIN_DIR)/%.o,$(CC_SRCS)) 


all: lib_tm.a Makefile

init:
					mkdir -p $(BIN_DIR)
					@echo $(C_OBJS)
					@echo $(CC_OBJS)

lib_tm.a: init $(C_OBJS) $(CC_OBJS)
					rm -f $@
					ar -qc $@ $(C_OBJS) $(CC_OBJS)

# For C source files
$(BIN_DIR)/%.o: $(SRC_DIR)/%.c Makefile
	${CXX} ${CXXFLAGS} -o $@ -c $<

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cc Makefile
	${CXX} ${CXXFLAGS} -o $@ -c $<

default: all 

array:		lib_tm.a Makefile
						$(MAKE) -C apps/array array

moldyn:		lib_tm.a Makefile
						$(MAKE) -C apps/moldyn moldyn

clean:
		-rm -f -R $(BIN_DIR)
		-rm -f lib_tm.a

clean_array:
		$(MAKE) -C apps/array clean

clean_moldyn:
		$(MAKE) -C apps/moldyn clean

clean_all:	clean clean_array clean_moldyn

.PHONY: all clean
