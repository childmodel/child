/***************************************************************************/
/**
**  @file
**  @brief unit test code for Tipper triangulator
**
** Use with:
**   (stand-alone)
**   g++ -DTIPPER_TEST -DDONT_USE_PREDICATE TipperTriangulator_test.cpp
**     TipperTriangulator.cpp TipperTriangulatorError.cpp
**   (with algorithms in predicates.cpp)
**   g++ -DTIPPER_TEST TipperTriangulator_test.cpp TipperTriangulator.cpp
**      TipperTriangulatorError.cpp ../Predicates/predicates.cpp
**
*/
/***************************************************************************/

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>

#include "TipperTriangulator.h"

#ifndef DONT_USE_PREDICATE
#include "../globalFns.h"
Predicates predicate;
#endif

// generate output files
const bool WRITE_FILES = false;

static
void sanity_check_ccwedge(int nedges, const edge* edges){
  // test ccw edge code
  for(int iedge=0; iedge<nedges; iedge++){
    const oriented_edge e1(iedge,true);
    const oriented_edge ccw_from = e1.ccw_edge_around_from(edges);
    if (ccw_from.o())
      assert( edges[iedge].from == edges[ccw_from.e()].from);
    else
      assert( edges[iedge].from == edges[ccw_from.e()].to);

    const oriented_edge e2(iedge,false);
    const oriented_edge ccw_to = e2.ccw_edge_around_from(edges);
    if (ccw_to.o())
      assert( edges[iedge].to == edges[ccw_to.e()].from);
    else
      assert( edges[iedge].to == edges[ccw_to.e()].to);

    if (0) //DEBUG
      cout << "edge=" << iedge
	   << " ret=" << edges[iedge].ret
	   << " lef=" << edges[iedge].lef
	   << " ccw_from=" << ccw_from.e()
	   << " ccw_to=" << ccw_to.e()
	   << endl;
  }
}

static
void sanity_check_edge(const edge *edges){
  // some sanity checks
  int i=0;
  while(edges[i].from != edge::end ) {
    int leftp = -1,rightp = -1;
    {
      const int& from = edges[i].from;
      const int& lef = edges[i].lef;
      const int& ref = edges[i].ref;
      if (lef != edge::none) {
	leftp = (edges[lef].from==from) ? edges[lef].to : edges[lef].from;
      }
      if (ref != edge::none) {
	rightp = (edges[ref].from==from) ? edges[ref].to : edges[ref].from;
      }
    }
    int leftp2 = -1,rightp2 = -1;
    {
      const int& to = edges[i].to;
      const int& let = edges[i].let;
      const int& ret = edges[i].ret;
      if (let != edge::none) {
	leftp2 = (edges[let].to==to) ? edges[let].from : edges[let].to;
      }
      if (ret != edge::none) {
	rightp2 = (edges[ret].to==to) ? edges[ret].from : edges[ret].to;
      }
    }
    assert(leftp==leftp2);
    if (rightp!=rightp2) {
      cout << "ERR edge=" << i
	   << " pto=" << edges[i].to << " pfrom=" << edges[i].from
	   << " left(from)=" << leftp << " left(to)=" << leftp2
	   << " righp(from)=" << rightp << " rightp(to)=" << rightp2
	   << endl;
    }
    assert(rightp==rightp2);

    i++;
  }
}

static
void write_point(int npoints, const point* p)
{
  ofstream file("points_sorted");
  file.precision(10);
  file << npoints << '\n';
  for(int i=0;i<npoints;++i){
    p[i].write(file);
  }
  file << flush;
}
static
void write_edge(const edge *edges, const point* p)
{
  ofstream file("triggy");
  file.precision(10);
  int i=0;
  while(edges[i].from != edge::end ){
    edges[i].write(file,p);
    file << "edge=" << i
	 << " pfrom=" << edges[i].from << " pto=" << edges[i].to
	 << " let=" << edges[i].let << " lef=" << edges[i].lef
	 << " ret=" << edges[i].ret << " ref=" << edges[i].ref
	 << '\n';
    i++;
  }
  file << flush;
}

static
void test_sort_triangulate(int npoints, point *p){
  int nedges;
  edge* edges = NULL;

  int npoints_unique;
  tt_sort_triangulate(npoints,p,&npoints_unique,&nedges,&edges);
  cout << "npoints=" << npoints
       << " npoints_unique=" << npoints_unique << endl;

  if (WRITE_FILES)
    write_point(npoints, p);
  if (WRITE_FILES)
    write_edge(edges, p);

  sanity_check_edge(edges);
  sanity_check_ccwedge(nedges, edges);

  int nelem;
  elem* elems = NULL;

  tt_build_elem_table(npoints_unique, p, nedges, edges, &nelem, &elems);
  cout << "nelem=" << nelem << endl;

  delete [] elems;

  delete [] edges;
}

static
void generate_dataset(int n,point* p){
  //make a set of points perturbed from a uniform grid
  srand(0);
  for (int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      p[i+j*n] = point(double(i)+double(rand())/RAND_MAX*n*1.e-3,
		       double(j)+double(rand())/RAND_MAX*n*1.e-3);
    }
  }
  if (WRITE_FILES) {
    const long npoints=n*n;
    ofstream file("points_res");
    file << npoints << '\n';
    for(int i=0;i<npoints;++i){
      p[i].write(file);
    }
    file << flush;
  }
}

static
void test_triangulate_random(int n){
  const long npoints=n*n;
  //set up the point array
  point *p = new point[npoints];
  generate_dataset(n,p);

  test_sort_triangulate(npoints,p);

  delete [] p;
}

void test_triangulate_from_file(){
  long npoints;
  ifstream file("points");
  assert( file.good() );
  file >> npoints;
  point *p = new point[npoints];
  for(int i=0;i<npoints;++i){
    double x, y, z;
    int b;
    file >> x;
    file >> y;
    file >> z; // ignored
    file >> b; // ignored
    p[i] = point(x,y,i+1);
  }

  test_sort_triangulate(npoints,p);

  delete [] p;
}

#if !defined(DONT_USE_MAIN)
int main(void){

  switch(2){
  case 0:
    {
      int n=3;
      cout << "n= " << n << endl;
      test_triangulate_random(n);
    }
    break;
  case 1:
    {
      for(int n=2; n!=100; n++){
	cout << "n= " << n << endl;
	test_triangulate_random(n);
      }
    }
    break;
  default:
    {
      test_triangulate_from_file();
    }
    break;
  }
  return 0;
}
#endif
