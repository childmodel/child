//A triangulation routine based on Tipper's convex hull algorithm
//Computers and geoscience vol17 no 5 pp 597-632,1991
//Scaling is nlogn for random datasets
//Mike Bithell 31/08/01


// test code for triangulation 

#include <math.h>
#include <fstream.h>
#include <stdlib.h>
#include <assert.h>

#include "hulltr_e.h"
#include "hulltr_test.h"

void sanity_check_edge(edge *edges){
  // some sanity check
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
    //cout << i
    //   << " " << leftp << " " << leftp2 << " " << rightp << " " << rightp2
    //   << endl;
    assert(leftp==leftp2);
    if (rightp!=rightp2) {
      cout << "ERR edge=" << i
	   << " pto=" << edges[i].to << " pfrom=" << edges[i].from
	   << " left(from)=" << leftp << " left(to)=" << leftp2
	   << " righp(from)=" << rightp << " rightp(to)=" << rightp2
	   << endl;
    } else {
#if 0
      cout << "edge=" << i
	   << " pto=" << edges[i].to << " pfrom=" << edges[i].from
	   << " left(from)=" << leftp << " left(to)=" << leftp2
	   << " righp(from)=" << rightp << " rightp(to)=" << rightp2
	   << endl;
#endif
    }
    assert(rightp==rightp2);

    i++;      
  }
}

void write_point(int npoints, point* p)
{
  ofstream file("points_sorted");
  file << npoints << endl;
  for(int i=0;i<npoints;++i){
    p[i].write(file);
  }
}
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


//#define WRITE_FILES

void test_sort_triangulate(int npoints, point *p){
  int nedges;
  edge* edges = NULL;

  sort_triangulate(npoints,p,&nedges,&edges);

#if defined(WRITE_FILES)
  write_point(npoints, p);
#endif
#if defined(WRITE_FILES)
  write_edge(edges, p);
#endif

#if 1
  sanity_check_edge(edges);
#endif

  int nelem;
  elem* elems = NULL;

  build_elem_table(npoints, p, nedges, edges, &nelem, &elems);
  cout << "nelem=" << nelem << endl;

  delete [] elems;

  delete [] edges;
}

void generate_dataset(int n,point* p){
  //make a set of points perturbed from a uniform grid
  srand(0);
  for (int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      p[i+j*n].x=double(i)+double(rand())/RAND_MAX*n*1.e-3;
      p[i+j*n].y=double(j)+double(rand())/RAND_MAX*n*1.e-3;
    } 
  }
#if defined(WRITE_FILES)
  {
    const long npoints=n*n;
    ofstream file("points");
    file << npoints << endl;
    for(int i=0;i<npoints;++i){
      p[i].write(file);
    }
  }
#endif
}

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

#if defined(USE_MAIN)
int main(){
  //  const int n=100; // makes it fail on Sun
#if 1
#if 0
  int n=3;
  cout << "n= " << n << endl;
  test_triangulate_random(n);
  
#else
  for(int n=2; n!=100; n++){
    cout << "n= " << n << endl;
    test_triangulate_random(n);
  }
#endif
#else
  test_triangulate_from_file();
#endif
}
#endif
