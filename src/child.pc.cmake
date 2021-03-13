Name: child
Description: Child library
Version: ${CHILD_VERSION}
Libs: -L${CMAKE_INSTALL_PREFIX}/lib -lchild
Cflags: -I${CMAKE_INSTALL_PREFIX}/include -I${CMAKE_INSTALL_PREFIX}/include/child -I${CMAKE_INSTALL_PREFIX}/include/child/ChildInterface

