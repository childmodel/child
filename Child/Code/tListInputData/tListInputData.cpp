/**************************************************************************/
/**
**  @file tListInputData.cpp
**  @brief Functions for class tListInputData.
**
**  Modifications:
**   - changed .tri file format from points-edges-triangles to
**     points-triangles-edges, compatible with earlier format (gt 1/98)
**   - GT merged tListIFStreams and tListInputData into a single class
**     to avoid multiple definition errors resulting from mixing
**     template & non-template classes (1/99)
**   - Bug fix in constructor: nnodes was being read from edge and
**     triangle files -- thus arrays dimensioned incorrectly! (GT 04/02)
**
**  $Id: tListInputData.cpp,v 1.12 2003-05-16 14:01:06 childcvs Exp $
*/
/**************************************************************************/

#include "tListInputData.h"
#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif

template< class tSubNode >
void tListInputData< tSubNode >::
ReportIOError(IOErrorType t, const char *filename,
	      const char *suffix, int n) const{
  cerr << "\nFile: '" << filename << suffix << "' "
       << "- Can't read ";
  switch (t){
  case IOTime:
    cerr << "time";
    break;
  case IOSize:
    cerr << "size";
    break;
  case IORecord:
    cerr << "record " << n;
    break;
  }
  cerr << "." << endl;
  ReportFatalError( "Input/Output Error." );
}


/**************************************************************************\
**
**  tListInputData constructor
**
**  The constructor does most of the work. It takes an input file (i.e.,
**  the "main" input file) and reads from it the base name of the files
**  that contain the triangulation. It then opens <basename>.nodes,
**  <basename>.z, <basename>.edges, and <basename>.tri. Assuming the
**  files are valid, the desired time-slice is read from infile, and
**  the start of data for that  time-slice is sought in each of the four
**  triangulation files. The arrays are dimensioned as needed, and
**  GetFileEntry() is called to read the data into the arrays. Note that
**  the time in each file is identified by a space character preceding it
**  on the same line.
**
**  Modifications:
**   - Calls to seekg appear not to be working properly with the g++
**     compiler, and a quick web search revealed all sorts of complaints
**     about seekg under g++. So I recoded the file reading to avoid
**     using seekg. (GT Feb 01)
**   - Fixed bug in which no. edges and triangles were incorrectly
**     assigned to nnodes, instead of nedges and ntri. (GT 04/02)
**
\**************************************************************************/
template< class tSubNode >
tListInputData< tSubNode >::
tListInputData( tInputFile &infile )                   //tListInputData
{
   int righttime;                   // flag: we've found the right time slice
   double time, intime;             // current & desired time
   char basename[80],               // base name of input files
       inname[80];                  // full name of an input file
   char headerLine[kMaxNameLength]; // header line read from input file

   // Read base name for triangulation files from infile
   infile.ReadItem( basename, "INPUTDATAFILE" );
   
   // Open each of the four files
   strcpy( inname, basename );
   strcat( inname, SNODES );
   nodeinfile.open(inname);    // Node input file pointer
   //assert( nodeinfile.good() );
   strcpy( inname, basename );
   strcat( inname, SEDGES );
   edgeinfile.open(inname);    // Edge input file pointer
   //assert( edgeinfile.good() );
   strcpy( inname, basename );
   strcat( inname, STRI );
   triinfile.open( inname );   // Triangle input file pointer
   //assert( triinfile.good() );
   strcpy( inname, basename );
   strcat( inname, SZ );
   zinfile.open( inname );     // Elevations input file pointer
   //assert( zinfile.good() );

   // Make sure we found them
   if( !nodeinfile.good() || !edgeinfile.good() || !triinfile.good()
       || !zinfile.good() )
   {
      cerr << "Error: I can't find one or more of the following files:\n"
           << "\t" << basename << SNODES "\n"
           << "\t" << basename << SEDGES "\n"
           << "\t" << basename << STRI "\n"
           << "\t" << basename << SZ "\n";
      ReportFatalError( "Unable to open triangulation input file(s)." );
   }
   
   // Find out which time slice we want to extract
   intime = infile.ReadItem( intime, "INPUTTIME" );
   cout << "intime = " << intime << endl;
   cout << "Is node input file ok? " << nodeinfile.good()
	<< " Are we at eof? " << nodeinfile.eof() << endl;

   // Find specified input times in input data files and read # items.
   // First, nodes:
   righttime = 0;
   while( !( nodeinfile.eof() ) && !righttime )
   {
     /*nodeinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         nodeinfile.seekg( -nodeinfile.gcount(), ios::cur );
         nodeinfile >> time;
         cout << "from file, time = " << time << endl;
	 cout << "Read: " << headerLine << endl;
         if( time == intime ) righttime = 1;
      }*/
     nodeinfile >> time;
     if (nodeinfile.fail())
       ReportIOError(IOTime, basename, SNODES);
     cout << "Read time: " << time << endl;
     if( time != intime )
       {
	 int i;
	 nodeinfile >> nnodes;
	 if (nodeinfile.fail())
	   ReportIOError(IOSize, basename, SNODES);
	 cout << "nnodes = " << nnodes << endl;
	 for( i=1; i<=nnodes+1; i++ ) {
	   nodeinfile.getline( headerLine, kMaxNameLength );
	 }
       }
     else righttime = 1;
      cout << " NOW are we at eof? " << nodeinfile.eof() << endl;
   }
   if( !( nodeinfile.eof() ) ) {
     nodeinfile >> nnodes;
     if (nodeinfile.fail())
       ReportIOError(IOSize, basename, SNODES);
   } else
   {
      cerr << "Couldn't find the specified input time in the node file\n";
      ReportFatalError( "Input error" );
   }
   // Then elevations (or "z" values):
   righttime = 0;
   while( !( zinfile.eof() ) && !righttime )
   {
     /*zinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         zinfile.seekg( -zinfile.gcount(), ios::cur );
         zinfile >> time;
         if( time == intime ) righttime = 1;
      }*/
     zinfile >> time;
     if (zinfile.fail())
       ReportIOError(IOTime, basename, SZ);
     cout << "Read time: " << time << endl;
     if( time != intime )
       {
	 int i;
	 zinfile >> nnodes;
	 if (zinfile.fail())
	   ReportIOError(IOSize, basename, SZ);
	 cout << "nnodes = " << nnodes << endl;
	 for( i=1; i<=nnodes+1; i++ ) {
	   zinfile.getline( headerLine, kMaxNameLength );
	 }
       }
     else righttime = 1;  
   }
   if( !( zinfile.eof() ) ) {
     zinfile >> nnodes;
     if (zinfile.fail())
       ReportIOError(IOSize, basename, SZ);
   } else
   {
      cerr << "Couldn't find specified input time in elevation file\n";
      ReportFatalError( "Input error" );
   }
   // Now edges:
   righttime = 0;
   while( !( edgeinfile.eof() ) && !righttime )
   {
     /*edgeinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         edgeinfile.seekg( -edgeinfile.gcount(), ios::cur );
         edgeinfile >> time;
         if( time == intime ) righttime = 1;
      }*/

     edgeinfile >> time;
     if (edgeinfile.fail())
       ReportIOError(IOTime, basename, SEDGES);
     cout << "Read time: " << time << endl;
     if( time != intime )
       {
	 int i;
	 edgeinfile >> nedges;
	 if (edgeinfile.fail())
	   ReportIOError(IOSize, basename, SEDGES);
	 cout << "nedges = " << nedges << endl;
	 for( i=1; i<=nedges+1; i++ ) {
	   edgeinfile.getline( headerLine, kMaxNameLength );
	 }
       }
     else righttime = 1;  
   }
   if( !( edgeinfile.eof() ) ) {
     edgeinfile >> nedges;
     if (edgeinfile.fail())
       ReportIOError(IOSize, basename, SEDGES);
   } else
   {
      cerr << "Couldn't find the specified input time in the edge file\n";
      ReportFatalError( "Input error" );
   }

   // And finally, triangles:
   righttime = 0;
   while( !( triinfile.eof() ) && !righttime )
   {
     /*triinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         triinfile.seekg( -triinfile.gcount(), ios::cur );
         triinfile >> time;
         if( time == intime ) righttime = 1;
      }*/
     triinfile >> time;
     if (triinfile.fail())
       ReportIOError(IOTime, basename, STRI);
     cout << "Read time: " << time << endl;
     if( time != intime )
       {
	 int i;
	 triinfile >> ntri;
	 if (triinfile.fail())
	   ReportIOError(IOSize, basename, STRI);
	 cout << "ntri = " << ntri << endl;
	 for( i=1; i<=ntri+1; i++ ) {
	   triinfile.getline( headerLine, kMaxNameLength );
	 }
       }
     else righttime = 1;  

   }
   if( !( triinfile.eof() ) ) {
     triinfile >> ntri;
     if (triinfile.fail())
       ReportIOError(IOSize, basename, STRI);
   } else
   {
      cerr << "Couldn't find the specified input time in the tri file\n";
      ReportFatalError( "Input error" );
   }

   // Dimension the arrays accordingly
   x.setSize( nnodes );
   y.setSize( nnodes );
   z.setSize( nnodes );
   edgid.setSize( nnodes );
   boundflag.setSize( nnodes );
   orgid.setSize( nedges );
   destid.setSize( nedges );
   nextid.setSize( nedges );
   p0.setSize( ntri );
   p1.setSize( ntri );
   p2.setSize( ntri );
   e0.setSize( ntri );
   e1.setSize( ntri );
   e2.setSize( ntri );
   t0.setSize( ntri );
   t1.setSize( ntri );
   t2.setSize( ntri );
   
   // Read in data from file
   GetFileEntry();
   
   // Close the files
   nodeinfile.close();
   edgeinfile.close();
   triinfile.close();
   zinfile.close();
   
}


/**************************************************************************\
**
**  tListInputData::GetFileEntry
**
**  Reads node, edge, and triangle data from the four triangulation input
**  files. Assumes that each files is open and valid and that the current
**  reading point in each corresponds the start of data for the desired
**  time-slice.
**
\**************************************************************************/
template< class tSubNode >
void tListInputData< tSubNode >::
GetFileEntry()                  //tListInputData
{
   char const * const  basename = "<file>";
   int i;

   for( i=0; i< nnodes; i++ ){
       nodeinfile >> x[i] >> y[i] >> edgid[i] >> boundflag[i];
       if (nodeinfile.fail())
	 ReportIOError(IORecord, basename, SNODES, i);
       zinfile >> z[i];
       if (zinfile.fail())
	 ReportIOError(IORecord, basename, SZ, i);
   }
   
   for( i=0; i<nedges; i++ ) {
       edgeinfile >> orgid[i] >> destid[i] >> nextid[i];
       if (edgeinfile.fail())
	 ReportIOError(IORecord, basename, SEDGES, i);
   }
   for( i=0; i< ntri; i++ ) {
       triinfile >> p0[i] >> p1[i] >> p2[i] >> t0[i] >> t1[i] >> t2[i] 
                 >> e0[i] >> e1[i] >> e2[i];
       if (triinfile.fail())
	 ReportIOError(IORecord, basename, STRI, i);
   }
}


/**************************************************************************\
**
**  tListInputData::GetKeyEntry
**
**  Provides alternative keyboard entry of triangulation data for 
**  testing purposes. Not currently supported.
**
\**************************************************************************/
template< class tSubNode >
void tListInputData< tSubNode >::
GetKeyEntry()                   //tListInputData
{
   int i;
   for( i=0; i < nnodes; i++ )
   {
      cout << "x y z edgid boundary:" << endl;
      cin >> x[i] >> y[i] >> z[i] >> edgid[i] >> boundflag[i];
   }
   for( i=0; i < nedges; i++ )
   {
      cout << "orgid destid nextid" << endl;
      cin >> orgid[i] >> destid[i] >> nextid[i];
   }
   for( i=0; i< ntri; i++ )
   {
      cout << "nodeids (3), edgids (3), triangleids (3)" << endl;
      cin >> p0[i] >> p1[i] >> p2[i]
          >> e0[i] >> e1[i] >> e2[i]
          >> t0[i] >> t1[i] >> t2[i];
   }

}


// code below this line is obsolete and can be deleted (gt 1/99)

/*
template< class tSubNode >
void tListInputData< tSubNode >::
GetFileEntry()                  //tListInputData
{
   int i;
   for( i=0; i< inchan.getNNodes(); i++ )
   {
      inchan.getNodeInFile() >> x[i] >> y[i]
                             >> edgid[i] >> boundflag[i];
      inchan.getZInFile() >> z[i];
   }
   
   for( i=0; i< inchan.getNEdges(); i++ )
       inchan.getEdgeInFile() >> orgid[i] >> destid[i] >> nextid[i];
   
   for( i=0; i< inchan.getNTri(); i++ )
       inchan.getTriInFile() >> p0[i] >> p1[i] >> p2[i]
                             >> t0[i] >> t1[i] >> t2[i]
                             >> e0[i] >> e1[i] >> e2[i];
   
}*/

/*
template< class tSubNode >
tListInputData< tSubNode >::
~tListInputData()                                                   //tInput
{
   cout << "    ~tListInputData()" << endl;
}*/

/**************************************************************************\
**
**         Utilities for classes tListIFStreams &  tListInputData
**
\**************************************************************************/
/*tListIFStreams::tListIFStreams()                               //tListIFStreams
{
   cout << "enter no. nodes" << endl;
   cin >> nnodes;
   cout << "enter no. edges" << endl;
   cin >> nedges;
   cout << "enter no. triangles" << endl;
   cin >> ntri;
   cout << "tInputChannel()" << endl;
}*/


/*
tListIFStreams::tListIFStreams( tInputFile &infile )
{
   int righttime;
   char thestring[80], inname[80];
   char headerLine[kMaxNameLength];
   infile.ReadItem( thestring, "INPUTDATAFILE" );
   
   strcpy( inname, thestring );
   strcat( inname, SNODES );
   nodeinfile.open(inname); // Node input file pointer
   assert( nodeinfile.good() );

   strcpy( inname, thestring );
   strcat( inname, SEDGES );
   edgeinfile.open(inname); // Edge input file pointer
   assert( edgeinfile.good() );

   strcpy( inname, thestring );
   strcat( inname, STRI );
   triinfile.open( inname ); //Triangle input file pointer
   assert( triinfile.good() );
   
   strcpy( inname, thestring );
   strcat( inname, SZ );
   zinfile.open( inname ); //Elevations input file pointer
   assert( zinfile.good() );
   
   intime = infile.ReadItem( intime, "INPUTTIME" );
   cout << "intime = " << intime << endl;
   //find specified input times in input data files and read no. items.
   //nodes:
   righttime = 0;
   while( !( nodeinfile.eof() ) && !righttime )
   {
      nodeinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         nodeinfile.seekg( -nodeinfile.gcount(), ios::cur );
         nodeinfile >> time;
         cout << "from file, time = " << time << endl;
         if( time == intime ) righttime = 1;
      }
   }
   if( !( nodeinfile.eof() ) ) nodeinfile >> nnodes;
   else
   {
      cerr << "Couldn't find the specified input time in the node file" << endl;
      ReportFatalError( "Input error" );
   }
   //elevations:
   righttime = 0;
   while( !( zinfile.eof() ) && !righttime )
   {
      zinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         zinfile.seekg( -zinfile.gcount(), ios::cur );
         zinfile >> time;
         if( time == intime ) righttime = 1;
      }
   }
   if( !( zinfile.eof() ) ) zinfile >> nnodes;
   else
   {
      cerr << "Couldn't find specified input time in elevation file" << endl;
      ReportFatalError( "Input error" );
   }
   //edges:
   righttime = 0;
   while( !( edgeinfile.eof() ) && !righttime )
   {
      edgeinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         edgeinfile.seekg( -edgeinfile.gcount(), ios::cur );
         edgeinfile >> time;
         if( time == intime ) righttime = 1;
      }
   }
   if( !( edgeinfile.eof() ) ) edgeinfile >> nedges;
   else
   {
      cerr << "Couldn't find the specified input time in the edge file" << endl;
      ReportFatalError( "Input error" );
   }
   //triangles:
   righttime = 0;
   while( !( triinfile.eof() ) && !righttime )
   {
      triinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         triinfile.seekg( -triinfile.gcount(), ios::cur );
         triinfile >> time;
         if( time == intime ) righttime = 1;
      }
   }
   if( !( triinfile.eof() ) ) triinfile >> ntri;
   else
   {
      cerr << "Couldn't find the specified input time in the tri file" << endl;
      ReportFatalError( "Input error" );
   }
}*/

/*
tListIFStreams::tListIFStreams( const char * argv = 0 )        //tListIFStreams
{
   if( argv == 0 )
   {
      cout << "enter no. nodes" << endl;
      cin >> nnodes;
      cout << "enter no. edges" << endl;
      cin >> nedges;
      cout << "enter no. triangles" << endl;
      cin >> ntri;
      //nodeinfile = 0;
   }
   else
   {
      char inname[80];
      strcpy( inname, argv );
      strcat( inname, SNODES );
      nodeinfile.open(inname); // Node input file pointer
      assert( nodeinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, SEDGES );
      edgeinfile.open(inname); // Edge input file pointer
      assert( edgeinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, STRI );
      triinfile.open( inname ); //Triangle input file pointer
      assert( triinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, SZ );
      zinfile.open( inname ); //Elevations input file pointer
      assert( zinfile.good() );
      
      nodeinfile >> time >> nnodes;
      edgeinfile >> time >> nedges;
      triinfile >> time >> ntri;
      zinfile >> time >> nnodes;
   }
   cout << "tListIFStreams( argv )" << endl;
}*/

/*
tListIFStreams::~tListIFStreams()                              //tListIFStreams
{   
   nodeinfile.close();
   edgeinfile.close();
   triinfile.close();
   zinfile.close();
   cout << "    ~tListIFStreams()" << endl;
}

int tListIFStreams::getNNodes() const {return nnodes;}        //tListIFStreams

int tListIFStreams::getNEdges() const {return nedges;}        //tListIFStreams

int tListIFStreams::getNTri() const {return ntri;}            //tListIFStreams

ifstream &tListIFStreams::getNodeInFile()                     //tListIFStreams
{return nodeinfile;}

ifstream &tListIFStreams::getEdgeInFile()                     //tListIFStreams
{return edgeinfile;}

ifstream &tListIFStreams::getTriInFile()                      //tListIFStreams
{return triinfile;}

ifstream &tListIFStreams::getZInFile()                      //tListIFStreams
{return zinfile;}
*/

//default constructor
//template< class tSubNode >
//tListInputData< tSubNode >::tListInputData()                //tListInputData
//constructor
/*template< class tSubNode >
tListInputData< tSubNode >::
tListInputData( tListIFStreams &chan ) 
        : inchan( chan ),
          x( inchan.getNNodes() ), y( inchan.getNNodes() ),
          z( inchan.getNNodes() ), edgid( inchan.getNNodes() ),
          boundflag( inchan.getNNodes() ),
          orgid( inchan.getNEdges() ), destid( inchan.getNEdges() ),
          nextid( inchan.getNEdges() ),
          p0( inchan.getNTri() ), p1( inchan.getNTri() ), p2( inchan.getNTri() ),
          e0( inchan.getNTri() ), e1( inchan.getNTri() ), e2( inchan.getNTri() ),
          t0( inchan.getNTri() ), t1( inchan.getNTri() ), t2( inchan.getNTri() )
{
   if( !inchan.getNodeInFile().good() ) GetKeyEntry( inchan );
   else GetFileEntry( inchan );
   cout << "tListInputData( inchan )" << endl;
}*/

/*
template< class tSubNode >
tListInputData< tSubNode >::
tListInputData( tInputFile &infile )                   //tListInputData
        : inchan( infile ),
          x( inchan.getNNodes() ), y( inchan.getNNodes() ),
          z( inchan.getNNodes() ), edgid( inchan.getNNodes() ),
          boundflag( inchan.getNNodes() ),
          orgid( inchan.getNEdges() ), destid( inchan.getNEdges() ),
          nextid( inchan.getNEdges() ),
          p0( inchan.getNTri() ), p1( inchan.getNTri() ), p2( inchan.getNTri() ),
          e0( inchan.getNTri() ), e1( inchan.getNTri() ), e2( inchan.getNTri() ),
          t0( inchan.getNTri() ), t1( inchan.getNTri() ), t2( inchan.getNTri() )
{
   if( !inchan.getNodeInFile().good() ) GetKeyEntry();
   else GetFileEntry();
   cout << "tListInputData( infile )" << endl;
}*/

