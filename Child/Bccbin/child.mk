PT = ../Code
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
EXENAME = child.exe

OBJ = obj

OBJECTS = childmain.$(OBJ) erosion.$(OBJ) meshElements.$(OBJ) mathutil.$(OBJ) \
 tInputFile.$(OBJ) tLNode.$(OBJ) tRunTimer.$(OBJ) tStreamMeander.$(OBJ) meander.$(OBJ) \
tStorm.$(OBJ) tStreamNet.$(OBJ) tUplift.$(OBJ) errors.$(OBJ) tFloodplain.$(OBJ) \
tEolian.$(OBJ) globalFns.$(OBJ) predicates.$(OBJ) tVegetation.$(OBJ) tListInputData.$(OBJ) \
ParamMesh_t.$(OBJ) TipperTriangulator.$(OBJ) TipperTriangulatorError.$(OBJ)

all : $(EXENAME)
.PHONY : all clean

$(EXENAME): $(OBJECTS)
	$(CXX) $(LDFLAGS) -e$@ $(OBJECTS) $(LIBS)

erosion.$(OBJ): $(PT)/Erosion/erosion.cpp
	$(CXX) $(CFLAGS) $(PT)/Erosion/erosion.cpp

meshElements.$(OBJ): $(PT)/MeshElements/meshElements.cpp
	$(CXX) $(CFLAGS) $(PT)/MeshElements/meshElements.cpp

mathutil.$(OBJ): $(PT)/Mathutil/mathutil.cpp
	$(CXX) $(CFLAGS) $(PT)/Mathutil/mathutil.cpp

tInputFile.$(OBJ): $(PT)/tInputFile/tInputFile.cpp
	$(CXX) $(CFLAGS) $(PT)/tInputFile/tInputFile.cpp

tLNode.$(OBJ): $(PT)/tLNode/tLNode.cpp
	$(CXX) $(CFLAGS) $(PT)/tLNode/tLNode.cpp

tRunTimer.$(OBJ): $(PT)/tRunTimer/tRunTimer.cpp
	$(CXX) $(CFLAGS) $(PT)/tRunTimer/tRunTimer.cpp

tListInputData.$(OBJ): $(PT)/tListInputData/tListInputData.cpp
	$(CXX) $(CFLAGS) $(PT)/tListInputData/tListInputData.cpp

tStorm.$(OBJ): $(PT)/tStorm/tStorm.h
	$(CXX) $(CFLAGS) $(PT)/tStorm/tStorm.cpp

tStreamNet.$(OBJ): $(PT)/tStreamNet/tStreamNet.cpp
	$(CXX) $(CFLAGS) $(PT)/tStreamNet/tStreamNet.cpp

tUplift.$(OBJ): $(PT)/tUplift/tUplift.cpp
	$(CXX) $(CFLAGS) $(PT)/tUplift/tUplift.cpp

errors.$(OBJ): $(PT)/errors/errors.cpp
	$(CXX) $(CFLAGS) $(PT)/errors/errors.cpp

tFloodplain.$(OBJ): $(PT)/tFloodplain/tFloodplain.cpp
	$(CXX) $(CFLAGS) $(PT)/tFloodplain/tFloodplain.cpp

tEolian.$(OBJ): $(PT)/tEolian/tEolian.cpp
	$(CXX) $(CFLAGS) $(PT)/tEolian/tEolian.cpp

ParamMesh_t.$(OBJ): $(PT)/tMesh/ParamMesh_t.cpp
	$(CXX) $(CFLAGS) $(PT)/tMesh/ParamMesh_t.cpp

TipperTriangulator.$(OBJ): $(PT)/tMesh/TipperTriangulator.cpp
	$(CXX) $(CFLAGS) $(PT)/tMesh/TipperTriangulator.cpp

TipperTriangulatorError.$(OBJ): $(PT)/tMesh/TipperTriangulatorError.cpp
	$(CXX) $(CFLAGS) $(PT)/tMesh/TipperTriangulatorError.cpp

globalFns.$(OBJ): $(PT)/globalFns.cpp
	$(CXX) $(CFLAGS) $(PT)/globalFns.cpp

predicates.$(OBJ): $(PT)/Predicates/predicates.cpp
	$(CXX) $(CFLAGS) $(PT)/Predicates/predicates.cpp

tVegetation.$(OBJ): $(PT)/tVegetation/tVegetation.cpp
	$(CXX) $(CFLAGS) $(PT)/tVegetation/tVegetation.cpp

tStreamMeander.$(OBJ): $(PT)/tStreamMeander/tStreamMeander.cpp
	$(CXX) $(CFLAGS) $(PT)/tStreamMeander/tStreamMeander.cpp

meander.$(OBJ): $(PT)/tStreamMeander/meander.cpp
	$(CXX) $(CFLAGS) $(PT)/tStreamMeander/meander.cpp

childmain.$(OBJ): $(PT)/childmain.cpp
	$(CXX) $(CFLAGS) $(PT)/childmain.cpp

clean::
	rm -f *.$(OBJ)
	rm -f $(EXENAME)

# dependencies: headers and template implementation files.
# use, for instance:
# find ${CHILDCODE} -name '*.h' | xargs grep -n -e include | grep '\.cpp'
HFILES = \
	$(PT)/compiler.h \
	$(PT)/Classes.h \
	$(PT)/Definitions.h \
	$(PT)/Erosion/erosion.h \
	$(PT)/Geometry/geometry.h \
	$(PT)/Inclusions.h \
	$(PT)/Mathutil/mathutil.h \
	$(PT)/MeshElements/meshElements.h \
	$(PT)/Predicates/predicates.h \
	$(PT)/errors/errors.h \
	$(PT)/globalFns.h \
	$(PT)/tArray/tArray.cpp \
	$(PT)/tArray/tArray.h \
	$(PT)/tEolian/tEolian.h \
	$(PT)/tFloodplain/tFloodplain.h \
	$(PT)/tInputFile/tInputFile.h \
	$(PT)/tLNode/tLNode.h \
	$(PT)/tList/tList.h \
	$(PT)/tListInputData/tListInputData.h \
	$(PT)/tMatrix/tMatrix.h \
	$(PT)/tMesh/ParamMesh_t.h \
	$(PT)/tMesh/TipperTriangulator.h \
	$(PT)/tMesh/heapsort.h \
	$(PT)/tMesh/tMesh.cpp \
	$(PT)/tMesh/tMesh.h \
	$(PT)/tMesh/tMesh2.cpp \
	$(PT)/tMeshList/tMeshList.h \
	$(PT)/tOutput/tOutput.cpp \
	$(PT)/tOutput/tOutput.h \
	$(PT)/tPtrList/tPtrList.h \
	$(PT)/tRunTimer/tRunTimer.h \
	$(PT)/tStorm/tStorm.h \
	$(PT)/tStreamMeander/tStreamMeander.h \
	$(PT)/tStreamMeander/meander.h \
	$(PT)/tStreamNet/tStreamNet.h \
	$(PT)/tUplift/tUplift.h \
	$(PT)/tVegetation/tVegetation.h \
	$(PT)/tStreamMeander/tStreamMeander.h \
	$(PT)/trapfpe.h

ParamMesh_t.$(OBJ): $(HFILES)
TipperTriangulator.$(OBJ) : $(HFILES)
TipperTriangulatorError.$(OBJ) : $(HFILES)
erosion.$(OBJ): $(HFILES)
errors.$(OBJ): $(HFILES)
globalFns.$(OBJ): $(HFILES)
mathutil.$(OBJ): $(HFILES)
meshElements.$(OBJ): $(HFILES)
predicates.$(OBJ): $(HFILES)
tEolian.$(OBJ): $(HFILES)
tFloodplain.$(OBJ): $(HFILES)
tInputFile.$(OBJ): $(HFILES)
tLNode.$(OBJ): $(HFILES)
tListInputData.$(OBJ): $(HFILES)
tRunTimer.$(OBJ): $(HFILES)
tStorm.$(OBJ) : $(HFILES)
tStreamNet.$(OBJ): $(HFILES)
tUplift.$(OBJ): $(HFILES)
tVegetation.$(OBJ): $(HFILES)
tStreamMeander.$(OBJ): $(HFILES)
meander.$(OBJ): $(HFILES)
childmain.$(OBJ) : $(HFILES)


