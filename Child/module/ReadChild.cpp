
// Arnaud Desitter

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream.h>

// for stat
#include <sys/types.h>
#include <sys/stat.h>

// include files from $EXPLORERHOME/include
#include<cx/UI.h>
#include<cx/Pyramid.h>
#include<cx/cxPyramid.api.h>
#include<cx/DataAccess.h>
#include<cx/DataOps.h>
#include<cx/PortAccess.h>

// Iris Explorer user function
#include <ReadChild.uf.h>

// local
#include <cxSmartPtr.h>
#include <readChildData.h>

// type of floating point for lattice
#define FLOAT float
#define CXFLOAT cx_prim_float

//#define DEBUG1

// extract a pointer to string hold by the "File Name" widget
static
const char *GetFileName(void){
  cxErrorCode ec = cx_err_none;
  int port = cxInputPortOpen("File Name");
  assert( -1 != port );
  cxParameter* Wdgt = static_cast<cxParameter*>(cxInputDataGet(port));
  assert(  cx_param_string == cxParameterTypeGet( Wdgt, &ec ) );
  
  return *( static_cast<char**>(cxParameterValueGet( Wdgt, &ec )) );
}

// compute and allocate a filename with the last suffix stripped
static
int IsFileNameValid(const char *filename){
  if (filename == NULL || 0 == strcmp(filename,"")) {
    cerr << "Filename is not set." << endl;
    return -1;
  }
  struct stat statbuf;
  if (0 != stat(filename,&statbuf)){
    cerr << "Filename is not readable." << endl;
    return -1;
  }
  if (! S_ISREG(statbuf.st_mode)){
    cerr << "Filename is invalid." << endl;
    return -1;
  }
  return 0;
}

// compute and allocate a filename with the last suffix stripped
static
int myBasename(const char *filename, char **basename){
  
  const size_t len_f = strlen(filename);
  // find location of the last '.' before a '/' if exist.
  int index_p = -1;
  for(size_t j = len_f; j != 0; j--){
    if (filename[j] == '/' ) break;
    if (filename[j] == '.' ) {
      index_p = j;
      break;
    }
  }
  if (index_p == -1 || index_p == 0){
    return -1;
  }
  
  char *b = new char[index_p+1];
  strncpy(b, filename, index_p);
  b[index_p] = '\0';
  *basename = b;
  return 0;
}

// test existence of files
static
int stat_file(const char* b){
  struct stat statbuf;
  if (0 != stat(b,&statbuf)){
    cerr << b << " is not readable." << endl;
    return -1;
  }
  if (! S_ISREG(statbuf.st_mode)){
    cerr << b << " is invalid." << endl;
    return -1;
  }
  return 0;
}

int testBasename(const char *basename){
  if (basename == NULL || 0 == strcmp(basename,"")) return -1;

  const size_t len_b = strlen(basename);
  char *b = new char[len_b+10];
  strcpy(b,basename);
  // test for .nodes

  assert(b[len_b] == '\0'); // paranoid !
  b[len_b] = '\0';
  strcat(b,".nodes");
  if (0 != stat_file(b))
    goto fail;
  b[len_b] = '\0';
  strcat(b,".tri");
  if (0 != stat_file(b))
    goto fail;
  b[len_b] = '\0';
  strcat(b,".z");
  if (0 != stat_file(b))
    goto fail;
  
  delete [] b;
  return 0;
 fail:
  delete [] b;
  return -1;
}

void ReadChild
(
 int nStep,
 int TypeVariable,
 cxPyramid    **Pyr2D,
 cxParameter  **currentTime
 ) {
#ifdef DEBUG1
  cout << "Entering main function" << endl; 
#endif
  //-----------------------------------------------------------------
  // Extract File Name
  char *basename(0);
  {
    char const * const filename = GetFileName();
    // early return if filename is not set. No message.
    if (filename == NULL || 0 == strcmp(filename,""))
      return;
#ifdef DEBUG1
    cout << "FileName=" << filename << endl;
#endif
    if (0 != IsFileNameValid(filename))
      return;
    if (0 != myBasename(filename, &basename)){
      cerr << "Filename does not have a suffix." << endl;
      return;
    }
#ifdef DEBUG1
    cout << "BaseName=" << basename << endl;
#endif
    if (0 != testBasename(basename))
      return;
  }
  //-----------------------------------------------------------------
  ReadChildData ChildData;
  try {
    if (!ChildData.LoadData(basename, nStep, TypeVariable)) {
      cerr << "Time step \"" << nStep << "\" does not exist." << endl;
      delete [] basename;
      return;
    }
  } catch (const BadFile &ex) {
    cerr << "Error when reading file \"" << ex.filename()
	 << "\": " << ex.ExceptionName()
	 << "." << endl;
    delete [] basename;
    return;
  } catch (...) {
    cerr << "A unkown exception occured during file processing." << endl;
    delete [] basename;
    return;
  }
  delete [] basename; basename = 0;
  const size_t npoints = ChildData.nnodes();
  const size_t nelem = ChildData.ntri();
  assert( npoints > 0 );
  assert( nelem > 0 );

#ifdef DEBUG1
  cout << "Data loaded." << endl;
#endif
  // set time
  *currentTime = cxParamDoubleNew(ChildData.currentTime());
  //-----------------------------------------------------------------
  // Pyramid building
  cxLattice_var baselat;                // base lattice
  cxConnection_var conn;                // elt2node connectivity
  cxPyramid_var Pyr2Dt;
  FLOAT *latcoord  = 0;                 // lattice coordinate buffer
  FLOAT *latdata  = 0;                  // lattice data buffer
  long *elem_  = 0;
  long *connec  = 0;
  cxErrorCode ec = cx_err_none;
  
  long nbn = npoints;   // nb of nodes
  const long nbe = nelem;     // nb of elements
  const long nbnve = 3*nelem; // nb of vertices
                              // 3 nodes per triangles
    
#ifdef DEBUG1
  cerr << "Nodes "    << nbn   << endl 
       << "Elements " << nbe   << endl
       << "Vertices " << nbnve << endl;
#endif
  //-----------------------------------------------------------------
  // Allocate Explorer Type

  // Allocate base lattice
  baselat= cxLatNew(1,            /* 1D lattice */
		    &nbn,         /* where dims[] is nbne */
		    2,            /* variables / node */
		    CXFLOAT,      /* type is float */
		    2,            /* number of coordinates / node */
		    // 3 for x,y,z
		    cx_coord_curvilinear /* curvilinear coordinates */
		    );
  if (cxDataAllocErrorGet()){
    cxModAlert("Error. Memory allocation failed for base lattice.");
    cxDataAllocErrorClear();
    return;
  }
  // Allocate pyramid
  Pyr2Dt = cxPyrNew(0);
  if (cxDataAllocErrorGet()){
    cxModAlert("Error. Memory allocation failed for Pyramid.");
    cxDataAllocErrorClear();
    return;
  }
  // Allocate connectivity
  conn = cxConnNew(nbe,nbnve);
  if (cxDataAllocErrorGet()){
    cxModAlert("Error. Memory allocation failed for connectivity");
    cxDataAllocErrorClear();
    return;
  }
  // Allocate and set Pyramid Dictionary
  long *iptr  = 0;
  {
    cxPyramidDictionary_var dict;
    dict = cxPyrDictDefault( 2 );         // triangle/quadrangle dict.
    if ( !dict ) {
      cxModAlert("Error while allocating dictionnary");
      cxDataAllocErrorClear();
      return;
    }
    cxPyramidDictionarySet(Pyr2Dt, dict, &ec );
    if ( cx_err_none != ec ){
      cxModAlert("Error while setting dictionnary");
      cxDataAllocErrorClear();
      return;
    } else {
      dict._retn();
    }
  }
  cxPyramidCompressionTypeSet(Pyr2Dt,cx_compress_unique, &ec);
  if ( cx_err_none != ec ){
    cxModAlert("Error during cxPyramidCompressionTypeSet.");
    cxDataAllocErrorClear();
    return;
  }
  iptr = (long *)cxPyramidCompressionIndexGet(Pyr2Dt, &ec);
  if ( cx_err_none != ec ){
    cxModAlert("Error during cxPyramidCompressionIndexGet.");
    cxDataAllocErrorClear();
    return;
  }
  iptr[0] = cx_pyramid_dict_triangle;
    
    //-----------------------------------------------------------------
    // Build and set base Lattice 
    /* Fetch pointers and fill it with coordinates and data */
  {
    void *_latdata, *_latcoord;
    ec=cxLatPtrGet(baselat, 
		   NULL,
		   &_latdata,
		   NULL,
		   &_latcoord
		   );
    latdata = static_cast<FLOAT *>(_latdata);
    latcoord = static_cast<FLOAT *>(_latcoord);
  }
  if ( cx_err_none != ec){
    cxModAlert("Error. base lattice corrupted.");
    cxDataAllocErrorClear();
    return;
  }

  {
    for(int inode=0;inode<nbn;inode++){
      latcoord[2*inode  ] = ChildData.x[inode];
      latcoord[2*inode+1] = ChildData.y[inode];
    }
  }
  {
    // Data
    for(int inode=0;inode<nbn;inode++){
      latdata[2*inode  ] = ChildData.z[inode];
      latdata[2*inode+1] = ChildData.data[inode];
    }
  }

  ec = cxPyrSet(Pyr2Dt,baselat);   /* set base lattice */
  if ( cx_err_none != ec){
    cxModAlert("Error during set base lattice.");
    cxDataAllocErrorClear();
    return;
  } else {
    baselat._retn();
  }

  //-----------------------------------------------------------------
  // Set Pyramid Layer and set base Lattice by recovering data from
  // ReadChildData object

  ec = cxPyrLayerSet(Pyr2Dt,1,NULL,NULL); /* nullify first connectivity level */
  if ( cx_err_none != ec){
    cxModAlert("Error during cxPyrLayerSet.");
    cxDataAllocErrorClear();
    return;
  }
    
  /* Fetch pointers and set up the connections list */
  ec = cxConnPtrGet(conn,NULL,NULL,&elem_,&connec);
  if ( cx_err_none != ec){
    cxModAlert("Error during cxConnPtrGet.");
    cxDataAllocErrorClear();
    return;
  }
  {
    // element to vertex connectivity table
    for(int ielem=0;ielem<nbe;ielem++){
      connec[3*ielem  ]= ChildData.p0[ielem];
      connec[3*ielem+1]= ChildData.p1[ielem];
      connec[3*ielem+2]= ChildData.p2[ielem];
    }
  }
  // build element to node table
  for(int ielem=0;ielem<(nbe+1);ielem++){
    elem_[ielem]= ielem * 3;  /* 3 nodes per elements */
  }
  
  ec = cxPyrLayerSet(Pyr2Dt,2,conn,NULL);
  if ( cx_err_none != ec){
    cxModAlert("Error during cxPyrLayerSet.");
    cxDataAllocErrorClear();
    return;
  } else {
    conn._retn();
  }
  
  //-----------------------------------------------------------------
  // Return the pyramid
  *Pyr2D = Pyr2Dt._retn();
  
  return;
}
