/**************************************************************************\
**
**  tListInputData.cpp
**
**  Functions for classes tListInputData and tListIFStreams.
**
**  Modifications:
**   - changed .tri file format from points-edges-triangles to
**     points-triangles-edges, compatible with earlier format (gt 1/98)
**
**  $Id: tListInputData.cpp,v 1.4 1998-02-12 01:45:19 stlancas Exp $
\**************************************************************************/

#include "tListInputData.h"

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
tListIFStreams::tListIFStreams( tInputFile &infile )
{
   int righttime;
   char thestring[80], inname[80];
   char headerLine[kMaxNameLength];
   infile.ReadItem( thestring, "INPUTDATAFILE" );
   
   strcpy( inname, thestring );
   strcat( inname, ".nodes" );
   nodeinfile.open(inname); /* Node input file pointer */
   assert( nodeinfile.good() );

   strcpy( inname, thestring );
   strcat( inname, ".edges" );
   edgeinfile.open(inname); /* Edge input file pointer */
   assert( edgeinfile.good() );

   strcpy( inname, thestring );
   strcat( inname, ".tri" );
   triinfile.open( inname ); //Triangle input file pointer
   assert( triinfile.good() );
   
   strcpy( inname, thestring );
   strcat( inname, ".z" );
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
}

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
      strcat( inname, ".nodes" );
      nodeinfile.open(inname); /* Node input file pointer */
      assert( nodeinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, ".edges" );
      edgeinfile.open(inname); /* Edge input file pointer */
      assert( edgeinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, ".tri" );
      triinfile.open( inname ); //Triangle input file pointer
      assert( triinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, ".z" );
      zinfile.open( inname ); //Elevations input file pointer
      assert( zinfile.good() );
      
      nodeinfile >> time >> nnodes;
      edgeinfile >> time >> nedges;
      triinfile >> time >> ntri;
      zinfile >> time >> nnodes;
   }
   cout << "tListIFStreams( argv )" << endl;
}

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

template< class tSubNode >
tListInputData< tSubNode >::
tListInputData( tInputFile &infile )                             //tListInputData
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
}

template< class tSubNode >
void tListInputData< tSubNode >::
GetKeyEntry()                   //tListInputData
{
   int i;
   for( i=0; i < inchan.getNNodes(); i++ )
   {
      cout << "x y z edgid boundary:" << endl;
      cin >> x[i] >> y[i] >> z[i] >> edgid[i] >> boundflag[i];
   }
   for( i=0; i < inchan.getNEdges(); i++ )
   {
      cout << "orgid destid nextid" << endl;
      cin >> orgid[i] >> destid[i] >> nextid[i];
   }
   for( i=0; i< inchan.getNTri(); i++ )
   {
      cout << "nodeids (3), edgids (3), triangleids (3)" << endl;
      cin >> p0[i] >> p1[i] >> p2[i]
          >> e0[i] >> e1[i] >> e2[i]
          >> t0[i] >> t1[i] >> t2[i];
   }
}

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
   
}

template< class tSubNode >
tListInputData< tSubNode >::
~tListInputData()                                                   //tInput
{
   cout << "    ~tListInputData()" << endl;
}

