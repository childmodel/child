# GNU compiler

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

# -march=i686: generates code for pentiumpro and later
# -march=pentium3: generates code for pentium III and later
ARCH := -march=pentium3
# optimise
CFLAGS = $(WARNINGFLAGS) -g -O2 $(ARCH) -c
LDFLAGS = $(WARNINGFLAGS) -g -O2 $(ARCH)
# no optimisation, build is faster
#CFLAGS = $(WARNINGFLAGS) -g $(ARCH) -c
#LDFLAGS = $(WARNINGFLAGS) -g $(ARCH)
LIBS =

LDFLAGS += -o $@

OBJEXT = o
EXEEXT =
