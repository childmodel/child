PT = $(CHILDCODE)
CC = CC
WARNINGFLAGS = 
CFLAGS = $(WARNINGFLAGS) -c -O
LDFLAGS = $(WARNINGFLAGS) -O
LIBS = -L/usr/local/lib -lm

OBJECTS = toddlermain.o erosion.o meshElements.o mathutil.o \
 tInputFile.o tLNode.o tRunTimer.o \
tPtrList.o tStorm.o tStreamNet.o tUplift.o errors.o tFloodplain.o \
tEolian.o globalFns.o predicates.o tVegetation.o

toddler: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o toddler $(LIBS)

erosion.o: $(PT)/Erosion/erosion.cpp $(PT)/Erosion/erosion.h
	$(CC) $(CFLAGS) $(PT)/Erosion/erosion.cpp

meshElements.o: $(PT)/MeshElements/meshElements.cpp \
 $(PT)/MeshElements/meshElements.h
	$(CC) $(CFLAGS) $(PT)/MeshElements/meshElements.cpp

mathutil.o: $(PT)/Mathutil/mathutil.cpp $(PT)/Mathutil/mathutil.h
	$(CC) $(CFLAGS) $(PT)/Mathutil/mathutil.cpp

tInputFile.o: $(PT)/tInputFile/tInputFile.cpp $(PT)/tInputFile/tInputFile.h
	$(CC) $(CFLAGS) $(PT)/tInputFile/tInputFile.cpp

tLNode.o: $(PT)/tLNode/tLNode.cpp $(PT)/tLNode/tLNode.h
	$(CC) $(CFLAGS) $(PT)/tLNode/tLNode.cpp

tPtrList.o: $(PT)/tPtrList/tPtrList.cpp $(PT)/tPtrList/tPtrList.h
	$(CC) $(CFLAGS) $(PT)/tPtrList/tPtrList.cpp

tRunTimer.o: $(PT)/tRunTimer/tRunTimer.cpp $(PT)/tRunTimer/tRunTimer.h
	$(CC) $(CFLAGS) $(PT)/tRunTimer/tRunTimer.cpp

tStorm.o: $(PT)/tStorm/tStorm.h $(PT)/tStorm/tStorm.cpp
	$(CC) $(CFLAGS) $(PT)/tStorm/tStorm.cpp

tStreamNet.o: $(PT)/tStreamNet/tStreamNet.cpp $(PT)/tStreamNet/tStreamNet.h
	$(CC) $(CFLAGS) $(PT)/tStreamNet/tStreamNet.cpp

tUplift.o: $(PT)/tUplift/tUplift.cpp $(PT)/tUplift/tUplift.h
	$(CC) $(CFLAGS) $(PT)/tUplift/tUplift.cpp

errors.o: $(PT)/errors/errors.cpp $(PT)/errors/errors.h
	$(CC) $(CFLAGS) $(PT)/errors/errors.cpp

tFloodplain.o: $(PT)/tFloodplain/tFloodplain.cpp \
$(PT)/tFloodplain/tFloodplain.h
	$(CC) $(CFLAGS) $(PT)/tFloodplain/tFloodplain.cpp

tEolian.o: $(PT)/tEolian/tEolian.cpp $(PT)/tEolian/tEolian.h
	$(CC) $(CFLAGS) $(PT)/tEolian/tEolian.cpp

ParamMesh_t.o: $(PT)/tMesh/ParamMesh_t.cpp $(PT)/tMesh/ParamMesh_t.h
	$(CC) $(CFLAGS) $(PT)/tMesh/ParamMesh_t.cpp

globalFns.o: $(PT)/globalFns.cpp $(PT)/globalFns.h
	$(CC) $(CFLAGS) $(PT)/globalFns.cpp

predicates.o: $(PT)/Predicates/predicates.cpp $(PT)/Predicates/predicates.h
	$(CC) $(CFLAGS) $(PT)/Predicates/predicates.cpp

tVegetation.o: $(PT)/tVegetation/tVegetation.cpp $(PT)/tVegetation/tVegetation.h
	$(CC) $(CFLAGS) $(PT)/tVegetation/tVegetation.cpp

toddlermain.o: $(PT)/toddlermain.cpp
	$(CC) $(CFLAGS) $(PT)/toddlermain.cpp
