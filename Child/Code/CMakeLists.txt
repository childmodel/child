cmake_minimum_required(VERSION 2.8)

ENABLE_TESTING()
ADD_TEST(
  bmi_model_child_test
  ${CMAKE_CURRENT_BINARY_DIR}/bmi_model_child_test test_input_files.txt
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ChildInterface/tests/test_input_files.txt.cmake
  test_input_files.txt
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/child.pc.cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/child.pc
)

set(CMAKE_MACOSX_RPATH 1)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/Erosion
  ${CMAKE_CURRENT_SOURCE_DIR}/MeshElements
  ${CMAKE_CURRENT_SOURCE_DIR}/Mathutil
  ${CMAKE_CURRENT_SOURCE_DIR}/tIDGenerator
  ${CMAKE_CURRENT_SOURCE_DIR}/tInputFile
  ${CMAKE_CURRENT_SOURCE_DIR}/tLNode
  ${CMAKE_CURRENT_SOURCE_DIR}/tListInputData
  ${CMAKE_CURRENT_SOURCE_DIR}/tOption
  ${CMAKE_CURRENT_SOURCE_DIR}/tRunTimer
  ${CMAKE_CURRENT_SOURCE_DIR}/tStorm
  ${CMAKE_CURRENT_SOURCE_DIR}/tStratGrid
  ${CMAKE_CURRENT_SOURCE_DIR}/tTimeSeries
  ${CMAKE_CURRENT_SOURCE_DIR}/tStreamNet
  ${CMAKE_CURRENT_SOURCE_DIR}/tUplift
  ${CMAKE_CURRENT_SOURCE_DIR}/errors
  ${CMAKE_CURRENT_SOURCE_DIR}/tFloodplain
  ${CMAKE_CURRENT_SOURCE_DIR}/tEolian
  ${CMAKE_CURRENT_SOURCE_DIR}/tMesh
  ${CMAKE_CURRENT_SOURCE_DIR}/tMesh
  ${CMAKE_CURRENT_SOURCE_DIR}/tMesh
  ${CMAKE_CURRENT_SOURCE_DIR}/Predicates
  ${CMAKE_CURRENT_SOURCE_DIR}/tVegetation
  ${CMAKE_CURRENT_SOURCE_DIR}/tStreamMeander
  ${CMAKE_CURRENT_SOURCE_DIR}/tStreamMeander
  ${CMAKE_CURRENT_SOURCE_DIR}/tWaterSedTracker
  ${CMAKE_CURRENT_SOURCE_DIR}/tLithologyManager
)

set (CHILD_LIB_SRCS
  ChildInterface/bmi_model.cpp
  ChildInterface/child.cpp
  Erosion/erosion.cpp
  MeshElements/meshElements.cpp
  Mathutil/mathutil.cpp
  tIDGenerator/tIDGenerator.cpp
  tInputFile/tInputFile.cpp
  tLNode/tLNode.cpp
  tListInputData/tListInputData.cpp
  tOption/tOption.cpp
  tRunTimer/tRunTimer.cpp
  tStorm/tStorm.cpp
  tStratGrid/tStratGrid.cpp
  tTimeSeries/tTimeSeries.cpp
  tStreamNet/tStreamNet.cpp
  tUplift/tUplift.cpp
  errors/errors.cpp
  tFloodplain/tFloodplain.cpp
  tEolian/tEolian.cpp
  tMesh/ParamMesh_t.cpp
  tMesh/TipperTriangulator.cpp
  tMesh/TipperTriangulatorError.cpp
  globalFns.cpp
  Predicates/predicates.cpp
  tVegetation/tVegetation.cpp
  tStreamMeander/tStreamMeander.cpp
  tStreamMeander/meander.cpp
  tWaterSedTracker/tWaterSedTracker.cpp
  tLithologyManager/tLithologyManager.cpp
)

add_library(child-shared SHARED ${CHILD_LIB_SRCS})
install(TARGETS child-shared DESTINATION lib)
set_target_properties(child-shared PROPERTIES OUTPUT_NAME "child")

add_library(child-static STATIC ${CHILD_LIB_SRCS})
set_target_properties(child-static PROPERTIES OUTPUT_NAME "child")
install(TARGETS child-static DESTINATION lib)

SET(CHILD_EXE_SRCS CHILDINTERFACE/CHILDDRIVER.CPP)

add_executable (child ${CHILD_EXE_SRCS})
target_link_libraries (child child-static)

install(FILES child.pc DESTINATION lib/pkgconfig)

install(TARGETS child DESTINATION bin)

add_executable(bmi_model_child_test ChildInterface/tests/bmi_model_child_test.cpp)
target_link_libraries(bmi_model_child_test child-shared)

install(
  FILES ChildInterface/bmi_model.h ChildInterface/child.h
  DESTINATION include/child/ChildInterface
)
install(
  FILES ChildInterface/childInterface.h
  DESTINATION include/child/ChildInterface)
install(
  FILES Erosion/erosion.h
  DESTINATION include/child/Erosion)
install(
  FILES Geometry/geometry.h
  DESTINATION include/child/Geometry)
install(
  FILES Mathutil/mathutil.h
  DESTINATION include/child/Mathutil)
install(
  FILES tIDGenerator/tIDGenerator.h
  DESTINATION include/child/tIDGenerator)
install(
  FILES MeshElements/meshElements.h
  DESTINATION include/child/MeshElements)
install(
  FILES Predicates/predicates.h
  DESTINATION include/child/Predicates)
install(
  FILES errors/errors.h
  DESTINATION include/child/errors)
install(
  FILES tArray/tArray.h tArray/tArray2.h tArray/tArray.cpp
  DESTINATION include/child/tArray)
install(
  FILES tEolian/tEolian.h
  DESTINATION include/child/tEolian)
install(
  FILES tFloodplain/tFloodplain.h
  DESTINATION include/child/tFloodplain)
install(
  FILES tInputFile/tInputFile.h
  DESTINATION include/child/tInputFile)
install(
  FILES tLNode/tLNode.h
  DESTINATION include/child/tLNode)
install(
  FILES tList/tList.h tList/tListFwd.h
  DESTINATION include/child/tList)
install(
  FILES tListInputData/tListInputData.h
  DESTINATION include/child/tListInputData)
install(
  FILES tLithologyManager/tLithologyManager.h
  DESTINATION include/child/tLithologyManager)
install(
  FILES tMatrix/tMatrix.h tMatrix/tMatrix.cpp
  DESTINATION include/child/tMatrix)
install(
  FILES
    tMesh/ParamMesh_t.h
    tMesh/TipperTriangulator.h
    tMesh/heapsort.h
    tMesh/tMesh.h
    tMesh/tMesh.cpp
    tMesh/tMesh2.cpp
  DESTINATION include/child/tMesh)
install(
  FILES tMeshList/tMeshList.h
  DESTINATION include/child/tMeshList)
install(
  FILES tOption/tOption.h
  DESTINATION include/child/tOption)
install(
  FILES tOutput/tOutput.h tOutput/tOutput.cpp
  DESTINATION include/child/tOutput)
install(
  FILES tPtrList/tPtrList.h
  DESTINATION include/child/tPtrList)
install(
  FILES tRunTimer/tRunTimer.h
  DESTINATION include/child/tRunTimer)
install(
  FILES tStorm/tStorm.h
  DESTINATION include/child/tStorm)
install(
  FILES tStratGrid/tStratGrid.h
  DESTINATION include/child/tStratGrid)
install(
  FILES tStreamMeander/meander.h tStreamMeander/tStreamMeander.h
  DESTINATION include/child/tStreamMeander)
install(
  FILES tStreamNet/tStreamNet.h
  DESTINATION include/child/tStreamNet)
install(
  FILES tTimeSeries/tTimeSeries.h
  DESTINATION include/child/tTimeSeries)
install(
  FILES tUplift/tUplift.h
  DESTINATION include/child/tUplift)
install(
  FILES tVegetation/tVegetation.h
  DESTINATION include/child/tVegetation)
install(
  FILES tWaterSedTracker/tWaterSedTracker.h
  DESTINATION include/child/tWaterSedTracker)
install(
  FILES
    Classes.h
    Definitions.h
    trapfpe.h
    Inclusions.h
    compiler.h
    globalFns.h
    Template_model.h
  DESTINATION include/child
)
