#platform
file = Makefile.unix
ifeq ($(OS),Windows_NT)
	file = Makefile.windows
endif

#include
include $(file)