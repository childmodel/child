PT = ../Code
LINT = C:\\lint\\lint-nt -iC:\\lint std.lnt
LINTFLAGS = child.lnt
MODFLAGS =  -u -zero -oo\($@\)

APPFLAGS =

OBJECTS = childmain.lob erosion.lob meshElements.lob mathutil.lob \
 tInputFile.lob tLNode.lob tRunTimer.lob tStreamMeander.lob meander.lob \
tStorm.lob tStreamNet.lob tUplift.lob errors.lob tFloodplain.lob \
tEolian.lob globalFns.lob predicates.lob tVegetation.lob \
ParamMesh_t.lob TipperTriangulator.lob TipperTriangulatorError.lob

.PHONY : all clean project module

all : project
module : $(OBJECTS)

project: $(OBJECTS)
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(OBJECTS)

erosion.lob: $(PT)/Erosion/erosion.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/Erosion/erosion.cpp

meshElements.lob: $(PT)/MeshElements/meshElements.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/MeshElements/meshElements.cpp

mathutil.lob: $(PT)/Mathutil/mathutil.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/Mathutil/mathutil.cpp

tInputFile.lob: $(PT)/tInputFile/tInputFile.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tInputFile/tInputFile.cpp

tLNode.lob: $(PT)/tLNode/tLNode.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tLNode/tLNode.cpp

tRunTimer.lob: $(PT)/tRunTimer/tRunTimer.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tRunTimer/tRunTimer.cpp

tStorm.lob: $(PT)/tStorm/tStorm.h
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tStorm/tStorm.cpp

tStreamNet.lob: $(PT)/tStreamNet/tStreamNet.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tStreamNet/tStreamNet.cpp

tUplift.lob: $(PT)/tUplift/tUplift.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tUplift/tUplift.cpp

errors.lob: $(PT)/errors/errors.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/errors/errors.cpp

tFloodplain.lob: $(PT)/tFloodplain/tFloodplain.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tFloodplain/tFloodplain.cpp

tEolian.lob: $(PT)/tEolian/tEolian.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tEolian/tEolian.cpp

ParamMesh_t.lob: $(PT)/tMesh/ParamMesh_t.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tMesh/ParamMesh_t.cpp

TipperTriangulator.lob: $(PT)/tMesh/TipperTriangulator.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tMesh/TipperTriangulator.cpp

TipperTriangulatorError.lob: $(PT)/tMesh/TipperTriangulatorError.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tMesh/TipperTriangulatorError.cpp

globalFns.lob: $(PT)/globalFns.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/globalFns.cpp

predicates.lob: $(PT)/Predicates/predicates.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/Predicates/predicates.cpp

tVegetation.lob: $(PT)/tVegetation/tVegetation.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tVegetation/tVegetation.cpp

tStreamMeander.lob: $(PT)/tStreamMeander/tStreamMeander.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tStreamMeander/tStreamMeander.cpp

meander.lob: $(PT)/tStreamMeander/meander.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/tStreamMeander/meander.cpp

childmain.lob: $(PT)/childmain.cpp
	$(LINT) $(LINTFLAGS) $(APPFLAGS) $(MODFLAGS) $(PT)/childmain.cpp

clean::
	rm -f *.lob

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
	$(PT)/tListInputData/tListInputData.cpp \
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

ParamMesh_t.lob: $(HFILES)
TipperTriangulator.lob : $(HFILES)
TipperTriangulatorError.lob : $(HFILES)
erosion.lob: $(HFILES)
errors.lob: $(HFILES)
globalFns.lob: $(HFILES)
mathutil.lob: $(HFILES)
meshElements.lob: $(HFILES)
predicates.lob: $(HFILES)
tEolian.lob: $(HFILES)
tFloodplain.lob: $(HFILES)
tInputFile.lob: $(HFILES)
tLNode.lob: $(HFILES)
tRunTimer.lob: $(HFILES)
tStorm.lob : $(HFILES)
tStreamNet.lob: $(HFILES)
tUplift.lob: $(HFILES)
tVegetation.lob: $(HFILES)
tStreamMeander.lob: $(HFILES)
meander.lob: $(HFILES)
childmain.lob : $(HFILES)


