# Borland compiler (Windows version)

CXX = bcc32
#
WARNINGFLAGS = -g255 -w-8027 -w-8008 -w-8066 -w-8070

# use Pentium pro instructions: "-6"
ARCH := -6
# optimise
CFLAGS = $(WARNINGFLAGS) $(ARCH) -c
LDFLAGS = $(WARNINGFLAGS) $(ARCH)
# no optimisation, build is faster
#CFLAGS = $(WARNINGFLAGS) -v $(ARCH) -c
#LDFLAGS = $(WARNINGFLAGS) -v $(ARCH)
LIBS =

LDFLAGS += -e$@

OBJEXT = obj
EXEEXT = .exe
