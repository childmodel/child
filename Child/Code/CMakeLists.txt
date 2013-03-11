cmake_minimum_required(VERSION 2.8)

ENABLE_TESTING ()
ADD_TEST (bmi_model_child_test ${CMAKE_CURRENT_BINARY_DIR}/bmi_model_child_test test_input_files.txt)
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ChildInterface/tests/test_input_files.txt.cmake test_input_files.txt)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/child.pc.cmake ${CMAKE_CURRENT_SOURCE_DIR}/child.pc )

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

set (child_LIB_SRCS
  ChildInterface/bmi_model_child.cpp
  ChildInterface/child.cpp
  ChildInterface/childDriver.cpp
  ChildInterface/childInterface.cpp
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

add_library (child-shared SHARED ${child_LIB_SRCS})
install (TARGETS child-shared DESTINATION lib COMPONENT child)
set_target_properties (child-shared PROPERTIES OUTPUT_NAME "child")

add_library (child-static STATIC ${child_LIB_SRCS})
set_target_properties (child-static PROPERTIES OUTPUT_NAME "child")
install (TARGETS child-static DESTINATION lib COMPONENT child)

set (child_SRCS ChildInterface/childDriver.cpp)

add_executable (child ${child_SRCS})
target_link_libraries (child child-static)

install(FILES child.pc DESTINATION lib/pkgconfig  COMPONENT child)

install (TARGETS child DESTINATION bin COMPONENT child)

add_executable (bmi_model_child_test ChildInterface/tests/bmi_model_child_test.cpp)
target_link_libraries (bmi_model_child_test child-shared)

install (FILES
  ChildInterface/bmi_model_child.h ChildInterface/child.h
  DESTINATION include/child/ChildInterface COMPONENT child)
install (FILES
  ChildInterface/childInterface.h
  DESTINATION include/child/ChildInterface COMPONENT child)
install (FILES
  Erosion/erosion.h
  DESTINATION include/child/Erosion COMPONENT child)
install (FILES
  Geometry/geometry.h
  DESTINATION include/child/Geometry COMPONENT child)
install (FILES
  Mathutil/mathutil.h
  DESTINATION include/child/Mathutil COMPONENT child)
install (FILES
  tIDGenerator/tIDGenerator.h
  DESTINATION include/child/tIDGenerator COMPONENT child)
install (FILES
  MeshElements/meshElements.h
  DESTINATION include/child/MeshElements COMPONENT child)
install (FILES
  Predicates/predicates.h
  DESTINATION include/child/Predicates COMPONENT child)
install (FILES
  errors/errors.h
  DESTINATION include/child/errors COMPONENT child)
install (FILES
  tArray/tArray.h
  tArray/tArray2.h
  tArray/tArray.cpp
  DESTINATION include/child/tArray COMPONENT child)
install (FILES
  tEolian/tEolian.h
  DESTINATION include/child/tEolian COMPONENT child)
install (FILES
  tFloodplain/tFloodplain.h
  DESTINATION include/child/tFloodplain COMPONENT child)
install (FILES
  tInputFile/tInputFile.h
  DESTINATION include/child/tInputFile COMPONENT child)
install (FILES
  tLNode/tLNode.h
  DESTINATION include/child/tLNode COMPONENT child)
install (FILES
  tList/tList.h
  tList/tListFwd.h
  DESTINATION include/child/tList COMPONENT child)
install (FILES
  tListInputData/tListInputData.h
  DESTINATION include/child/tListInputData COMPONENT child)
install (FILES
  tLithologyManager/tLithologyManager.h
  DESTINATION include/child/tLithologyManager COMPONENT child)
install (FILES
  tMatrix/tMatrix.h
  tMatrix/tMatrix.cpp
  DESTINATION include/child/tMatrix COMPONENT child)
install (FILES
  tMesh/ParamMesh_t.h
  tMesh/TipperTriangulator.h
  tMesh/heapsort.h
  tMesh/tMesh.h
  tMesh/tMesh.cpp
  tMesh/tMesh2.cpp
  DESTINATION include/child/tMesh COMPONENT child)
install (FILES
  tMeshList/tMeshList.h
  DESTINATION include/child/tMeshList COMPONENT child)
install (FILES
  tOption/tOption.h
  DESTINATION include/child/tOption COMPONENT child)
install (FILES
  tOutput/tOutput.h
  tOutput/tOutput.cpp
  DESTINATION include/child/tOutput COMPONENT child)
install (FILES
  tPtrList/tPtrList.h
  DESTINATION include/child/tPtrList COMPONENT child)
install (FILES
  tRunTimer/tRunTimer.h
  DESTINATION include/child/tRunTimer COMPONENT child)
install (FILES
  tStorm/tStorm.h
  DESTINATION include/child/tStorm COMPONENT child)
install (FILES
  tStratGrid/tStratGrid.h
  DESTINATION include/child/tStratGrid COMPONENT child)
install (FILES
  tStreamMeander/meander.h
  tStreamMeander/tStreamMeander.h
  DESTINATION include/child/tStreamMeander COMPONENT child)
install (FILES
  tStreamNet/tStreamNet.h
  DESTINATION include/child/tStreamNet COMPONENT child)
install (FILES
  tTimeSeries/tTimeSeries.h
  DESTINATION include/child/tTimeSeries COMPONENT child)
install (FILES
  tUplift/tUplift.h
  DESTINATION include/child/tUplift COMPONENT child)
install (FILES
  tVegetation/tVegetation.h
  DESTINATION include/child/tVegetation COMPONENT child)
install (FILES
  tWaterSedTracker/tWaterSedTracker.h
  DESTINATION include/child/tWaterSedTracker COMPONENT child)
install (FILES
  Classes.h
  Definitions.h
  trapfpe.h
  Inclusions.h
  compiler.h
  globalFns.h
  Template_model.h
  DESTINATION include/child COMPONENT child
)


