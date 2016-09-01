PT = ../Code

#include gccmac.mk
include gcc.mk
#include icc.mk
#include bcc.mk
#include PC-lint/pclint.mk

EXENAME = child$(EXEEXT)

OBJECTS = \
 childDriver.$(OBJEXT) \
 childInterface.$(OBJEXT) erosion.$(OBJEXT) \
 meshElements.$(OBJEXT) mathutil.$(OBJEXT) tIDGenerator.$(OBJEXT) \
 tInputFile.$(OBJEXT) tLNode.$(OBJEXT) tRunTimer.$(OBJEXT) \
 tStreamMeander.$(OBJEXT) meander.$(OBJEXT) \
 tStorm.$(OBJEXT) tStreamNet.$(OBJEXT) tUplift.$(OBJEXT) errors.$(OBJEXT) \
 tFloodplain.$(OBJEXT) tEolian.$(OBJEXT) globalFns.$(OBJEXT) \
 predicates.$(OBJEXT) tVegetation.$(OBJEXT) tListInputData.$(OBJEXT) \
 tStratGrid.$(OBJEXT) tOption.$(OBJEXT) \
 tTimeSeries.$(OBJEXT) ParamMesh_t.$(OBJEXT) TipperTriangulator.$(OBJEXT) \
 TipperTriangulatorError.$(OBJEXT) tWaterSedTracker.$(OBJEXT) \
 tLithologyManager.$(OBJEXT)

all : $(EXENAME)
.PHONY : all clean

$(EXENAME): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LIBS)

erosion.$(OBJEXT): $(PT)/Erosion/erosion.cpp
	$(CXX) $(CFLAGS) $(PT)/Erosion/erosion.cpp

meshElements.$(OBJEXT): $(PT)/MeshElements/meshElements.cpp
	$(CXX) $(CFLAGS) $(PT)/MeshElements/meshElements.cpp

mathutil.$(OBJEXT): $(PT)/Mathutil/mathutil.cpp
	$(CXX) $(CFLAGS) $(PT)/Mathutil/mathutil.cpp

tIDGenerator.$(OBJEXT): $(PT)/tIDGenerator/tIDGenerator.cpp
	$(CXX) $(CFLAGS) $(PT)/tIDGenerator/tIDGenerator.cpp

tInputFile.$(OBJEXT): $(PT)/tInputFile/tInputFile.cpp
	$(CXX) $(CFLAGS) $(PT)/tInputFile/tInputFile.cpp

tLNode.$(OBJEXT): $(PT)/tLNode/tLNode.cpp
	$(CXX) $(CFLAGS) $(PT)/tLNode/tLNode.cpp

tListInputData.$(OBJEXT): $(PT)/tListInputData/tListInputData.cpp
	$(CXX) $(CFLAGS) $(PT)/tListInputData/tListInputData.cpp

tLithologyManager.$(OBJEXT): $(PT)/tLithologyManager/tLithologyManager.cpp
	$(CXX) $(CFLAGS) $(PT)/tLithologyManager/tLithologyManager.cpp

tOption.$(OBJEXT): $(PT)/tOption/tOption.cpp
	$(CXX) $(CFLAGS) $(PT)/tOption/tOption.cpp

tRunTimer.$(OBJEXT): $(PT)/tRunTimer/tRunTimer.cpp
	$(CXX) $(CFLAGS) $(PT)/tRunTimer/tRunTimer.cpp

tStorm.$(OBJEXT): $(PT)/tStorm/tStorm.cpp
	$(CXX) $(CFLAGS) $(PT)/tStorm/tStorm.cpp

tStratGrid.$(OBJEXT):  $(PT)/tStratGrid/tStratGrid.cpp
	$(CXX) $(CFLAGS) $(PT)/tStratGrid/tStratGrid.cpp

tTimeSeries.$(OBJEXT): $(PT)/tTimeSeries/tTimeSeries.cpp
	$(CXX) $(CFLAGS) $(PT)/tTimeSeries/tTimeSeries.cpp

tStreamNet.$(OBJEXT): $(PT)/tStreamNet/tStreamNet.cpp
	$(CXX) $(CFLAGS) $(PT)/tStreamNet/tStreamNet.cpp

tUplift.$(OBJEXT): $(PT)/tUplift/tUplift.cpp
	$(CXX) $(CFLAGS) $(PT)/tUplift/tUplift.cpp

errors.$(OBJEXT): $(PT)/errors/errors.cpp
	$(CXX) $(CFLAGS) $(PT)/errors/errors.cpp

tFloodplain.$(OBJEXT): $(PT)/tFloodplain/tFloodplain.cpp
	$(CXX) $(CFLAGS) $(PT)/tFloodplain/tFloodplain.cpp

tEolian.$(OBJEXT): $(PT)/tEolian/tEolian.cpp
	$(CXX) $(CFLAGS) $(PT)/tEolian/tEolian.cpp

ParamMesh_t.$(OBJEXT): $(PT)/tMesh/ParamMesh_t.cpp
	$(CXX) $(CFLAGS) $(PT)/tMesh/ParamMesh_t.cpp

TipperTriangulator.$(OBJEXT): $(PT)/tMesh/TipperTriangulator.cpp
	$(CXX) $(CFLAGS) $(PT)/tMesh/TipperTriangulator.cpp

TipperTriangulatorError.$(OBJEXT): $(PT)/tMesh/TipperTriangulatorError.cpp
	$(CXX) $(CFLAGS) $(PT)/tMesh/TipperTriangulatorError.cpp

globalFns.$(OBJEXT): $(PT)/globalFns.cpp
	$(CXX) $(CFLAGS) $(PT)/globalFns.cpp

predicates.$(OBJEXT): $(PT)/Predicates/predicates.cpp
	$(CXX) $(CFLAGS) $(PT)/Predicates/predicates.cpp

tVegetation.$(OBJEXT): $(PT)/tVegetation/tVegetation.cpp
	$(CXX) $(CFLAGS) $(PT)/tVegetation/tVegetation.cpp

tWaterSedTracker.$(OBJEXT): $(PT)/tWaterSedTracker/tWaterSedTracker.cpp
	$(CXX) $(CFLAGS) $(PT)/tWaterSedTracker/tWaterSedTracker.cpp

tStreamMeander.$(OBJEXT): $(PT)/tStreamMeander/tStreamMeander.cpp
	$(CXX) $(CFLAGS) $(PT)/tStreamMeander/tStreamMeander.cpp

meander.$(OBJEXT): $(PT)/tStreamMeander/meander.cpp
	$(CXX) $(CFLAGS) $(PT)/tStreamMeander/meander.cpp

childInterface.$(OBJEXT): $(PT)/ChildInterface/childInterface.cpp
	$(CXX) $(CFLAGS) $(PT)/ChildInterface/childInterface.cpp

childDriver.$(OBJEXT): $(PT)/ChildInterface/childDriver.cpp
	$(CXX) $(CFLAGS) $(PT)/ChildInterface/childDriver.cpp

clean::
	rm -f $(EXENAME)
	rm -f *.$(OBJEXT)

# dependencies: headers and template implementation files.
# use, for instance:
# find ${CHILDCODE} -name '*.h' | xargs grep -n -e include | grep '\.cpp'
HFILES = \
	$(PT)/ChildInterface/childInterface.h \
	$(PT)/Classes.h \
	$(PT)/Definitions.h \
	$(PT)/Erosion/erosion.h \
	$(PT)/Geometry/geometry.h \
	$(PT)/Inclusions.h \
	$(PT)/Mathutil/mathutil.h \
	$(PT)/tIDGenerator/tIDGenerator.h \
	$(PT)/MeshElements/meshElements.h \
	$(PT)/Predicates/predicates.h \
	$(PT)/compiler.h \
	$(PT)/errors/errors.h \
	$(PT)/globalFns.h \
	$(PT)/tArray/tArray.cpp \
	$(PT)/tArray/tArray.h \
	$(PT)/tArray/tArray2.h \
	$(PT)/tEolian/tEolian.h \
	$(PT)/tFloodplain/tFloodplain.h \
	$(PT)/tInputFile/tInputFile.h \
	$(PT)/tLNode/tLNode.h \
	$(PT)/tList/tList.h \
	$(PT)/tList/tListFwd.h \
	$(PT)/tListInputData/tListInputData.h \
	$(PT)/tLithologyManager/tLithologyManager.h \
	$(PT)/tMatrix/tMatrix.h \
	$(PT)/tMesh/ParamMesh_t.h \
	$(PT)/tMesh/TipperTriangulator.h \
	$(PT)/tMesh/heapsort.h \
	$(PT)/tMesh/tMesh.cpp \
	$(PT)/tMesh/tMesh.h \
	$(PT)/tMesh/tMesh2.cpp \
	$(PT)/tMeshList/tMeshList.h \
	$(PT)/tOption/tOption.h \
	$(PT)/tOutput/tOutput.cpp \
	$(PT)/tOutput/tOutput.h \
	$(PT)/tPtrList/tPtrList.h \
	$(PT)/tRunTimer/tRunTimer.h \
	$(PT)/tStorm/tStorm.h \
	$(PT)/tStratGrid/tStratGrid.h \
	$(PT)/tStreamMeander/meander.h \
	$(PT)/tStreamMeander/tStreamMeander.h \
	$(PT)/tStreamNet/tStreamNet.h \
	$(PT)/tTimeSeries/tTimeSeries.h \
	$(PT)/tUplift/tUplift.h \
	$(PT)/tVegetation/tVegetation.h \
	$(PT)/tWaterSedTracker/tWaterSedTracker.h \
	$(PT)/trapfpe.h

ParamMesh_t.$(OBJEXT): $(HFILES)
TipperTriangulator.$(OBJEXT) : $(HFILES)
TipperTriangulatorError.$(OBJEXT) : $(HFILES)
childInterface.$(OBJEXT) : $(HFILES)
childInterfaceDriver.$(OBJEXT) : $(HFILES)
erosion.$(OBJEXT): $(HFILES)
errors.$(OBJEXT): $(HFILES)
globalFns.$(OBJEXT): $(HFILES)
mathutil.$(OBJEXT): $(HFILES)
meander.$(OBJEXT): $(HFILES)
meshElements.$(OBJEXT): $(HFILES)
predicates.$(OBJEXT): $(HFILES)
tEolian.$(OBJEXT): $(HFILES)
tFloodplain.$(OBJEXT): $(HFILES)
tInputFile.$(OBJEXT): $(HFILES)
tLNode.$(OBJEXT): $(HFILES)
tListInputData.$(OBJEXT): $(HFILES)
tLithologyManager.$(OBJEXT): $(HFILES)
tOption.$(OBJEXT): $(HFILES)
tRunTimer.$(OBJEXT): $(HFILES)
tStorm.$(OBJEXT) : $(HFILES)
tStratGrid.$(OBJEXT) : $(HFILES)
tStreamMeander.$(OBJEXT): $(HFILES)
tStreamNet.$(OBJEXT): $(HFILES)
tTimeSeries.$(OBJEXT) : $(HFILES)
tUplift.$(OBJEXT): $(HFILES)
tVegetation.$(OBJEXT): $(HFILES)
tWaterSedTracker.$(OBJEXT): $(HFILES)
