# Intel compiler

CXX = icpc
#
# Intel compiler is picky !
# 279: controlling expression is constant
# 981: operands are evaluated in unspecified order
# 810: conversion may loose significant bits
# 444: destructor for base class is not virtual
# 383: value copied to temporary, reference to temporary used
WARNINGFLAGS = -Wall -wd279,981,810 -wd444 -wd383
#WARNINGFLAGS =

# -march=pentiumiii -xK: generates code for pentium III and later
ARCH := -march=pentiumiii -xK
# optimise
PROFILE =# -p
CFLAGS = $(WARNINGFLAGS) -g $(PROFILE) -O2 $(ARCH) -c
LDFLAGS = $(WARNINGFLAGS) -g $(PROFILE) -O2 $(ARCH)
# no optimisation, build is faster
#CFLAGS = $(WARNINGFLAGS) -g -c
#LDFLAGS = $(WARNINGFLAGS) -g
LIBS =

LDFLAGS += -o $@

OBJEXT = o
EXEEXT =
