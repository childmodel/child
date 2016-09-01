#include <readChildData.h>

#include <fstream.h>
#include <string.h>
#include <assert.h>

// dumb string class to avoid leaks due to exceptions
// should use std::string if only the compiler wasn't so old.
class MyString {
  char *ptr_;
  MyString(const MyString&);
  MyString& operator=(const MyString&);
  void destroy();
public:
  MyString(char *);
  explicit MyString(size_t);
  ~MyString();
  char *retn();
  char *ptr() {
    return ptr_;
  }
  const char *c_str() const {
    return ptr_;
  }
};
void MyString::destroy() {
  delete [] ptr_;
  ptr_ = 0;
}
MyString::MyString(char *p){
  ptr_ = p;
}
MyString::MyString(size_t l){
  ptr_ = new char[l+1];
  ptr_[l] = '\0';
}
MyString::~MyString(){
  destroy();
}
char *MyString::retn(){
  char *p = ptr_;
  ptr_ = 0;
  return p;
}

enum {
  v_z = 0,
  v_slp,
  v_area,
  v_varea,
  v_q,
  v_qs,
  v_tau,
  v_tx
};

const char * const SuffixVariable[] =
{ 
  ".z",
  ".slp",
  ".area",
  ".varea",
  ".q",
  ".qs",
  ".tau",
  ".tx"
};

const size_t NSuffixVariable = sizeof(SuffixVariable)/sizeof(SuffixVariable[0]);

char const * const BadFileName[] = 
{
  "Can't open file",
  "File corrupted"
};

BadFile::BadFile(BadFileState_t s, const char *f) :
  filename_(0), state(s) 
{
  if (f) {
    filename_ = new char[strlen(f)+1];
    strcpy(filename_,f);
  }
}

BadFile::BadFile(const BadFile& b) :
  filename_(0), state(b.state) 
{
  if (b.filename_) {
    filename_ = new char[strlen(b.filename_)+1];
    strcpy(filename_,b.filename_);
  }
}

BadFile::~BadFile() {
  delete [] filename_;
}

const char *BadFile::filename() const {
  return filename_?filename_:"(unknown)";
}

const char *BadFile::ExceptionName() const {
  return BadFileName[state];
}

ReadChildData::ReadChildData() :
  nnodes_(0),
  ntri_(0),
  currentTime_(0.),
  x(0), y(0),
  p0(0), p1(0), p2(0),
  z(0), data(0)
{}

ReadChildData::~ReadChildData(){
  delete [] x;
  delete [] y;
  delete [] p0;
  delete [] p1;
  delete [] p2;
  delete [] z;
  delete [] data;
}

void ReadChildData::AllocateNodes(){
  delete [] x;
  x = new float[nnodes_];
  delete [] y;
  y = new float[nnodes_];
}

void ReadChildData::AllocateTriangles(){
  delete [] p0;
  p0 = new int[ntri_];
  delete [] p1;
  p1 = new int[ntri_];
  delete [] p2;
  p2 = new int[ntri_];
}

void ReadChildData::AllocateData(bool inZ){
  if (inZ) {
    delete [] z;
    z = new float[nnodes_];
  } else {
    delete [] data;
    data = new float[nnodes_];
  }
}

bool ReadChildData::LoadData(const char* basename, const int nStep,
			     const int TypeVariable)
{
  // preconditions
  assert(TypeVariable>=0 && TypeVariable < (int)NSuffixVariable);
  assert(basename != NULL);
  if (nStep < 1) {
    return false;
  }

  MyString filename(strlen(basename)+20);
  // node
  strcpy(filename.ptr(),basename);
  strcat(filename.ptr(),".nodes");
  if (!ReadNodes(filename.c_str(), nStep))
    goto fail;

  // element
  strcpy(filename.ptr(),basename);
  strcat(filename.ptr(),".tri");
  if (!ReadTriangles(filename.c_str(), nStep))
    goto fail;

  // elevation
  {
    const int zVariable = v_z;
    strcpy(filename.ptr(),basename);
    strcat(filename.ptr(),SuffixVariable[zVariable]);
    if (!ReadData(filename.c_str(), nStep, zVariable, true))
      goto fail;
  }

  // data
  if (TypeVariable == v_z){
    CopyZinData();
  } else {
    strcpy(filename.ptr(),basename);
    strcat(filename.ptr(),SuffixVariable[TypeVariable]);
    if (!ReadData(filename.c_str(), nStep, TypeVariable, false))
      goto fail;
  }

  return true;
 fail:
  return false;
}

void skipLine(ifstream& infile){
  char c;
  while(infile.get(c) && c != '\n');
}

// skip record until the "nStep"-th Step if reached
static
int skipRecords(ifstream& file, const int nStep, const char* filename){
  bool righttime = false;
  int StepRead = 0;
  while (! file.eof() && ! righttime) {
    StepRead++;
    if (StepRead < nStep) {
      skipLine(file);
      if ( file.bad() )
	throw BadFile(FileBad,filename);
      if (file.eof())
	continue;

      int n;
      file >> n;
      assert( !file.fail() );
      skipLine(file); // discards the rest of that line
      for(int i=0;i!=n;++i){
	skipLine(file);
	if ( file.bad() )
	  throw BadFile(FileBad,filename);
	if (file.eof())
	  continue;
      } 
    } else {
      righttime = true;
    }
  }
  if (righttime) 
    return 0;
  return -1;
}

void ReadChildData::CopyZinData(){
  AllocateData(false);
  memcpy(data,z,nnodes()*sizeof(data[0]));
}

bool ReadChildData::ReadNodes(const char* filename, const int nStep){
  ifstream nodeinfile;
  nodeinfile.open(filename);
  if (!nodeinfile.good()) {
    throw BadFile(CantOpen, filename);
  }

  if (0 != skipRecords(nodeinfile, nStep, filename)) {
    goto fail;
  }
  assert(nodeinfile.good());
  nodeinfile >> currentTime_;
  if ( nodeinfile.eof())
    goto fail;
  if (! nodeinfile.good())
    throw BadFile(FileBad,filename);
  size_t lnnodes;
  nodeinfile >> lnnodes;
  if ( ! nodeinfile.good() )
    throw BadFile(FileBad,filename);
  assert( !nodeinfile.fail() );
  nnodes_ = lnnodes;
  AllocateNodes();
  for(size_t i=0; i != nnodes_; ++i){
    int itemp1, itemp2;
    if (nodeinfile.eof()) throw BadFile(FileBad,filename);
    nodeinfile >> x[i] >> y[i] >> itemp1 >> itemp2;
    assert( !nodeinfile.fail() );
  }

  nodeinfile.close();
  return true;
 fail:
  nodeinfile.close();
  return false;
}

bool ReadChildData::ReadTriangles(const char* filename, const int nStep){
  ifstream triinfile;

  triinfile.open(filename);

  if (!triinfile.good()) {
    throw BadFile(CantOpen,filename);
  }

  if (0 != skipRecords(triinfile, nStep, filename)) {
    goto fail;
  }
  triinfile >> currentTime_;
  if ( triinfile.eof())
    goto fail;
  if (! triinfile.good())
    throw BadFile(FileBad,filename);
  triinfile >> ntri_;
  if ( ! triinfile.good() )
    throw BadFile(FileBad,filename);
  assert( !triinfile.fail() );
  AllocateTriangles();
  for(size_t i=0; i != ntri_; ++i){
    int
      itemp1, itemp2,itemp3,
      itemp4, itemp5,itemp6;
    if (triinfile.eof()) throw BadFile(FileBad,filename);

    triinfile >>
      p0[i] >> p1[i] >> p2[i] >>
      itemp1 >> itemp2  >> itemp3 >>
      itemp4 >> itemp5  >> itemp6;
    assert( !triinfile.fail() );
  }
  triinfile.close();
  return true;
 fail:
  triinfile.close();
  return false;
}

bool ReadChildData::ReadData(const char* filename, const int nStep,
			     int TypeVariable, bool inZ){
  // area is stored only for interior nodes
  const bool isArea = TypeVariable == v_area;
  ifstream datainfile;

  datainfile.open(filename);
  if (!datainfile.good()) {
    throw BadFile(CantOpen, filename);
  }

  if (0 != skipRecords(datainfile, nStep, filename)) {
    goto fail;
  }
  datainfile >> currentTime_;
  if ( datainfile.eof())
    goto fail;
  if (! datainfile.good())
    throw BadFile(FileBad,filename);

  size_t ntemp;
  datainfile >> ntemp;
  if ( ! datainfile.good() )
    throw BadFile(FileBad,filename);
  assert( !datainfile.fail() );
  if (isArea)
    assert(ntemp < nnodes_);
  else
    assert(ntemp == nnodes_);
  AllocateData(inZ);
  for(size_t i=0; i != ntemp; ++i){
    if (datainfile.eof()) throw BadFile(FileBad,filename);

    if (inZ)
      datainfile >> z[i];
    else
      datainfile >> data[i];
    assert( !datainfile.fail() );
  }
  if (isArea) {
    assert(!inZ);
    for(size_t i=ntemp; i != nnodes_; ++i){
      z[i] = 0.;
    }
  }
  datainfile.close();
  return true;
 fail:
  datainfile.close();
  return false;
}
