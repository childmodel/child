// test code for Tipper triangulator

#include <math.h>
#include <fstream.h>
#include <stdlib.h>
#include <assert.h>

#include "TipperTriangulator.h"

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
void sanity_check_edge(edge *edges){
  // some sanity checks
  int i=0;
  while(edges[i].from != -1 ) {
    int leftp = -1,rightp = -1;
    {
      const int from = edges[i].from;
      const int lef = edges[i].lef;
      const int ref = edges[i].ref;
      if (lef != -1) {
	if (edges[lef].from==from) leftp=edges[lef].to; else leftp=edges[lef].from;
      }
      if (ref != -1) {
	if (edges[ref].from==from) rightp=edges[ref].to; else rightp=edges[ref].from;
      }
    }
    int leftp2 = -1,rightp2 = -1;
    {
      const int to = edges[i].to;
      const int let = edges[i].let;
      const int ret = edges[i].ret;
      if (let != -1) {
	if (edges[let].to==to) leftp2=edges[let].from; else leftp2=edges[let].to;
      }
      if (ret != -1) {
	if (edges[ret].to==to) rightp2=edges[ret].from; else rightp2=edges[ret].to;
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
void write_point(int npoints, point* p)
{
  ofstream file("points_sorted");
  file << npoints << endl;
  for(int i=0;i<npoints;++i){
    p[i].write(file);
  }
}
static
void write_edge(edge *edges, point* p)
{
  ofstream file("triggy");
  int i=0;
  while(edges[i].from != -1 ){
    edges[i].write(file,p);
    file << "edge=" << i 
	 << " pfrom=" << edges[i].from << " pto=" << edges[i].to
	 << " let=" << edges[i].let << " lef=" << edges[i].lef 
	 << " ret=" << edges[i].ret << " ref=" << edges[i].ref 
	 << endl; 
    i++;
  }
}

static
void test_sort_triangulate(int npoints, point *p){
  int nedges;
  edge* edges = NULL;

  tt_sort_triangulate(npoints,p,&nedges,&edges);

  if (WRITE_FILES)
    write_point(npoints, p);
  if (WRITE_FILES)
    write_edge(edges, p);

  sanity_check_edge(edges);
  sanity_check_ccwedge(nedges, edges);

  int nelem;
  elem* elems = NULL;

  tt_build_elem_table(npoints, p, nedges, edges, &nelem, &elems);
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
      p[i+j*n].x=double(i)+double(rand())/RAND_MAX*n*1.e-3;
      p[i+j*n].y=double(j)+double(rand())/RAND_MAX*n*1.e-3;
    } 
  }
  if (WRITE_FILES) {
    const long npoints=n*n;
    ofstream file("points_res");
    file << npoints << endl;
    for(int i=0;i<npoints;++i){
      p[i].write(file);
    }
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
  file >> npoints;
  point *p = new point[npoints];
  for(int i=0;i<npoints;++i){
    file >> p[i].x;
    file >> p[i].y;
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
}
#endif
