
LINT = C:\\lint\\lint-nt -iC:\\lint std.lnt
LINTFLAGS = PC-lint\child.lnt
MODFLAGS =  -u -zero -oo\($@\)

APPFLAGS =

# map to compiler variable
CXX = $(LINT)
CFLAGS = $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS)
LDFLAGS = $(LINTFLAGS) $(APPFLAGS)
LIBS =

OBJEXT = lob

