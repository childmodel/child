# GNU compiler (mac OS X)

CXX = g++
#
# -O is necessary for -Wuninitialized to be on
# -Weffc++ -fmessage-length=0 gives useful but noisy warnings
WARNINGFLAGS = -pedantic -Wall -W \
	-Wwrite-strings \
	-Wpointer-arith -Wcast-qual -Wcast-align

# recent versions of Linux use new style casts in their headers
WARNINGFLAGS += -Wold-style-cast

# gcc 2.x does not put the standard C++ headers in the namespace "std"
# In such a case, remove the comment of the following line.
#WARNINGFLAGS += -DHAVE_NO_NAMESPACE

# -mpowerpc64 etc: generates code for 64 bit G5 (power pc 970)
# will not run on G3, G4, etc.
#ARCH := -mcpu=G5 -mtune=G5 -mpowerpc64 -mpowerpc-gpopt
ARCH :=
# optimise
#CFLAGS = $(WARNINGFLAGS) -g -O2 $(ARCH) -c
#LDFLAGS = $(WARNINGFLAGS) -g -O2 $(ARCH)
# no optimisation, build is faster
CFLAGS = $(WARNINGFLAGS) -g $(ARCH) -c
LDFLAGS = $(WARNINGFLAGS) -g $(ARCH)
LIBS =

LDFLAGS += -o $@

OBJEXT = o
EXEEXT =
