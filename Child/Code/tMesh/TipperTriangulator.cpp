//A triangulation routine based on Tipper's convex hull algorithm
//Computers and geoscience vol17 no 5 pp 597-632,1991
//Scaling is nlogn for random datasets
//Mike Bithell 31/08/01

// Arnaud Desitter. Q2 2002. Binding to CHILD data structures.

#include <math.h>
#include <fstream.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#define DEBUG_PRINT 1
#include "hulltr_e.h"

#if __SUNPRO_CC==0x420
typedef enum { false=0, true } bool;
#endif

#define TIMING 1

const point &point::operator=( const point &p ) {
  if ( &p != this ) {
    x=p.x; y=p.y;
  }
  return *this;
}
#if defined(DEBUG_PRINT)
void point::print () const {cout << x << ' '<< y <<endl;}
#endif
void point::write(ofstream& f) const {f<<x<<' '<<y<<endl;}


#if defined(DEBUG_PRINT)
void edge::print(const point p[]) const {p[from].print();p[to].print();}
#endif
void edge::write(ofstream& f,const point p[]) const {p[from].write(f);p[to].write(f);}
bool edge::visible(const point p[],int i) const {
  //test whether an edge on the hull is visible from a point
  //rely on the fact that a) hull is anticlockwise oriented
  //b)data is positive x ordered

  const double mindistance = 0.0000001;
  if (fabs(p[from].x-p[to].x) < mindistance) return true;
  if (fabs(p[to].y-p[from].y)<mindistance){
    if(p[from].x<p[to].x && p[i].y<p[from].y)return true;
    if(p[from].x>p[to].x && p[i].y>p[from].y)return true;
  }
  if (p[to].y>=p[i].y && p[from].y<=p[i].y && fabs(p[from].y-p[to].y)>mindistance)return true;

  if (p[to].x>p[from].x){
    if(p[i].y<p[from].y+(p[to].y-p[from].y)/(p[to].x-p[from].x)*(p[i].x-p[from].x))return true;
  }else if (p[to].x<p[from].x) {
    if(p[i].y>p[from].y+(p[to].y-p[from].y)/(p[to].x-p[from].x)*(p[i].x-p[from].x))return true;
  }
  return false;
}

int edge::swap(int tint,edge e[],const point p[]){
  
  //edge swapping routine - each edge has four neighbour edges
  //left and attached to from node (lef) right attached to to node (ret) etc.
  //these edges may be oriented so that their from node is that of 
  // the current edge, or not
  //this routine takes the lions share of the CPU - hull construction by comparison
  //takes much less (by about a factor of ten!)
  if (ref==-1 || lef==-1 || let==-1 || ret==-1)return 0;
  //test orientation of left and right edges - store the indices of
  //points that are not part of the current edge
  int leftp,rightp;
  if (e[lef].from==from)leftp=e[lef].to; else leftp=e[lef].from;
  if (e[ref].from==from)rightp=e[ref].to; else rightp=e[ref].from;
  const point 
    p1(p[leftp]-p[from]), p2(p[leftp]-p[to]), 
    p3(p[rightp]-p[from]), p4(p[rightp]-p[to]);
  double dt1=p1.dot(p2);double dt2=p3.dot(p4);
  //only do the square roots if we really need to - saves a bit of time
  if (dt1<0 || dt2<0){
    dt1=dt1/sqrt(p1.dot(p1)*p2.dot(p2));
    dt2=dt2/sqrt(p3.dot(p3)*p4.dot(p4));
    if ((dt1+dt2)<0){
      //now swap the left and right edges of neighbouring edges
      //taking into account orientation
      if (e[ref].from == from){
	e[ref].lef=lef;
	e[ref].let=tint;
      }else{
	e[ref].ref=tint;
	e[ref].ret=lef;
      }
      if (e[lef].from==from){
	e[lef].ref=ref;
	e[lef].ret=tint;
      }else{
	e[lef].lef=tint;
	e[lef].let=ref;
      }
      if (e[ret].to==to){
	e[ret].lef=tint;
	e[ret].let=let;
      }else{
	e[ret].ref=let;
	e[ret].ret=tint;
      }
      if (e[let].to==to){
	e[let].ref=tint;
	e[let].ret=ret;
      }else{
	e[let].lef=ret;
	e[let].let=tint;
      }
      //change the end-points for the current edge
      to=rightp;
      from=leftp;
      //re-jig the edges
      int rf=ref;
      ref=lef;
      int rt=ret;
      ret=rf;
      int lt=let;
      let=rt;
      lef=lt;
      //examine the neighbouring edges for delauniness recursively - this is
      //a lot more efficient than trying to swap all edges right at the end.
      e[lef].swap(lef,e,p);
      e[let].swap(let,e,p);
      e[ref].swap(ref,e,p);
      e[ret].swap(ret,e,p);
      return 1;
    }
  }
  return 0;
}

class cyclist{
  //a fixed size linked cyclical list using arrays
  //the code here is not yet robust necessarily to list shrinking to zero
  //number of elements
  const cyclist &operator=( const cyclist & );  // assignment operator
  cyclist( const cyclist & );
  void delNext(int list_pos);
  void add(int ej);
  void invariant() const;
  void checkRange(int e) const;
public:
  cyclist(int s);
  ~cyclist();
  int getEdge(int list_pos) const {
    checkRange(list_pos);
    return ejs[list_pos].data;
  }
  int delNextPos(int list_pos);
  int delNextNeg(int list_pos);
  int getNextPos(int list_pos) const {
    checkRange(list_pos);
    return ejs[list_pos].next;
  }
  int getNextNeg(int list_pos) const{
    checkRange(list_pos);
    return ejs[list_pos].prev;
  }
  int addBefore(int a, int ej);
  int addAfter(int a, int ej);
#if defined(DEBUG_PRINT)
  void print() const;
#endif
private:
  const int size; 
  struct item{
    int next,prev,data;
  };
  item* ejs;  // array of size "size"
  // hole points to the next unfilled location,
  // prev to the location that was last filled
  // we keep track of empty bits of the list using hole to point to an empty slot,
  // and the value ejs[hole] to point to the next empty slot
  int prev, hole, num;
};

cyclist::cyclist(int s): size(s),ejs(0),prev(0),hole(0),num(0) {
  ejs=new item[size];
  //fill ej array with a set of pointers to the next unfilled location
  for (int i=0;i<size;i++) ejs[i].data=i+1;
}
cyclist::~cyclist(){
  invariant();
  delete [] ejs;
}
void cyclist::checkRange(int e) const {
  assert(e<size);
  assert(0<=e);
}
void cyclist::invariant() const {
  assert(num<=size);
  assert(0<=num);
  checkRange(hole);
  checkRange(prev);
}
void cyclist::delNext(int list_pos){
  checkRange(list_pos);
  assert(num!=0);
  ejs[list_pos].data=hole;
  hole=list_pos;
  ejs[ejs[list_pos].prev].next=ejs[list_pos].next;
  ejs[ejs[list_pos].next].prev=ejs[list_pos].prev;
  num--;
}
int cyclist::delNextPos(int list_pos){
  delNext(list_pos);
  return ejs[list_pos].next;
}
int cyclist::delNextNeg(int list_pos){
  delNext(list_pos);
  return ejs[list_pos].prev;
}
void cyclist::add(int ej){
  //build hull from scratch in numerical order - we assume you got the orientation
  //right!! (anti-clockwise)
  checkRange(hole<size);
  const int n=ejs[hole].data;
  ejs[hole].data=ej;
  //prev stores the location of the place in the array that was most recently filled
  ejs[prev].next=hole;
  //rev is the set of backward pointers
  ejs[hole].prev=prev;
  //the list is cyclic
  ejs[hole].next=0;
  ejs[0].prev=hole;
  prev=hole;
  hole=n;num++;
}
int cyclist::addBefore(int a, int ej){
  checkRange(a);
  //first check for the empty list
  if (num ==0){add(ej);return prev;}
  //otherwise add on before the specified position, using the empty storage slot
  const int n=ejs[hole].data;
  ejs[hole].prev=ejs[a].prev;
  ejs[hole].next=a;
  ejs[ejs[a].prev].next=hole;
  ejs[a].prev=hole;
  ejs[hole].data=ej;
  prev=hole;
  hole=n;
  num++;
  //return the value that hole had at the start of the method 
  return prev;
}
int cyclist::addAfter(int a,int ej){
  checkRange(a);
  //first check for the empty list
  if (num ==0){add(ej);return prev;}
  //otherwise add on after the specified position, using the empty storage slot
  const int n=ejs[hole].data;
  ejs[hole].next=ejs[a].next;
  ejs[ejs[a].next].prev=hole;
  ejs[a].next=hole;
  ejs[hole].prev=a;
  ejs[hole].data=ej;
  prev=hole;
  hole=n;
  num++;
  //return the value that hole had at the start of the method 
  return prev;
}
#if defined(DEBUG_PRINT)
void cyclist::print() const {int j=ejs[0].next;for (int i=0;i<num;i++){cout<<ejs[j].data<<endl;j=ejs[j].next;}}
#endif

// vector product (or cross product) of p0p1,p0p2
static
double vecprod(int p0,int p1,int p2,const point *p){
  return
    (p[p1].x-p[p0].x)*(p[p2].y-p[p0].y)
    -(p[p1].y-p[p0].y)*(p[p2].x-p[p0].x);
}

void triangulate(int npoints,const point p[], int *pnedges, edge** edges_ret){

  //convex hull is a cyclical list - it will consist of anticlockwise
  //ordered edges - since each new point adds at most 1 extra edge (nett)
  //to the hull, there are at most npoints edges on the hull
  cyclist hull(npoints);
  
  //and the edges - there are at most three edges per point
  const long nn=3*npoints;
  edge* edges;
  edges=new edge[nn];

  //make first three edges  - these will form the initial convex hull
  //make sure orientation is anticlockwise
  // Arnaud: orientation is counter clockwise <=> angle(p0p1,p0p2) >= 0
  //  <=> sin(p0p1,p0p2) >=0 <=> vect_prod(p0p1,p0p2) >= 0
  {
    const double v = vecprod(0,1,2,p);
    if (v>0) {
      edges[0].from=0;
      edges[0].to=1;
      edges[1].from=1;
      edges[1].to=2;
      edges[2].from=2;
      edges[2].to=0;
    }
    else{
      edges[0].from=0;
      edges[0].to=2;
      edges[1].from=2;
      edges[1].to=1;
      edges[2].from=1;
      edges[2].to=0;
    }
  }
  //make left edges
  edges[0].lef=2;
  edges[0].let=1;
  edges[1].lef=0;
  edges[1].let=2;
  edges[2].lef=1;
  edges[2].let=0;
  // add the edges to the hull in order
  // get upper and lower edges as indices into the hull
  // 'cos hull is cyclic, we don't need to bother if these are visible
  // but they do need to be oriented so that upper edge is further round the 
  // hull in a positive direction than lower_hull_pos
  int start=hull.addAfter(0,0);
  int lower_hull_pos=hull.addAfter(start,1); 
  int upper_hull_pos=hull.addAfter(lower_hull_pos,2);
  // loop through the remaining points adding edges to the hull
  int next_edge=3;
  int count,saved_edge;
  
  for (int i=3;i<npoints;i++){
    saved_edge=-1;
    //go round the hull looking for visible edges - we need to go round
    //in two directions from the current upper and lower edges
    //first set up the new edge that joins to the point coincident
    //between upper and lower edge
    if(edges[hull.getEdge(upper_hull_pos)].visible(p,i)){
      //make new edge - from and to nodes preserve hull orientation
      edges[next_edge].from=edges[hull.getEdge(upper_hull_pos)].from;
      edges[next_edge].to=i;
      //save hull position for possible later use
      saved_edge=next_edge;
      //connectivity for swapping - we know the ID of the
      //next edge to be created since upper hull is visible
      edges[next_edge].lef=hull.getEdge(upper_hull_pos);
      edges[next_edge].let=next_edge+1;
      next_edge++;
    }
    else{
      //we can't see the upper edge - can we see the lower one?
      if(!edges[hull.getEdge(lower_hull_pos)].visible(p,i)){
	//can't see the upper or lower edge - chose a bad initial state! go round
	//the hull a bit.Upper hull will still be invisible (its what used to be lower hull)
	if (i==3){
	  lower_hull_pos=hull.getNextPos(upper_hull_pos);
	  upper_hull_pos=hull.getNextPos(lower_hull_pos);
	}else{
	  //or else its an error!
	  cout<<"Triangulate: Can't see the hull from the new point!? number is "<<i<<endl;
	  exit(1);
	}
      }
      edges[next_edge].from=i;
      edges[next_edge].to=edges[hull.getEdge(lower_hull_pos)].to;

      //connectivity for swapping - we know the ID of the
      //next edge but one to be created since upper hull is *not* visible
      edges[next_edge].let=hull.getEdge(lower_hull_pos);
      edges[next_edge].lef=next_edge+1;
      next_edge++;
    }
    //now we need to add the upper hull edges - drop through to set edge made above
    //to the new upper edge if upper edge not visible
    count=0;
    int h;
    while (edges[hull.getEdge(upper_hull_pos)].visible(p,i)){
      h=hull.getEdge(upper_hull_pos);
      edges[next_edge].from=i;
      edges[next_edge].to=edges[h].to;
      //if we got here, the upper hull is visible, so we know which way the edges point
      //set the left and right neighbour edges
      if (count!=0){
	//if we're on the second time round the loop, the edge made in
	// the previous pass needs to be connected on its right side
	edges[next_edge-1].ref=next_edge;
	edges[next_edge-1].ret=h;
      }
      count++;
      edges[h].ref=next_edge-1;
      edges[h].ret=next_edge;
      edges[next_edge].let=h;
      edges[next_edge].lef=next_edge-1;
      //check the hull edge's delauniness
      edges[h].swap(h,edges,p);
      next_edge++;
      //delete upper edge from the hull
      //and go round the hull in the positive direction
      upper_hull_pos=hull.delNextPos(upper_hull_pos);
    }
    //we drop through to here if upper edge is not visible
    //add edge to the hull
    //set the new upper edge to the most recently created upper edge
    upper_hull_pos=hull.addBefore(upper_hull_pos,next_edge-1);
    //now we need to add the lower hull edges
    while (edges[hull.getEdge(lower_hull_pos)].visible(p,i)){
      h=hull.getEdge(lower_hull_pos);
      //upper hull was not visible
      if (saved_edge==-1){
	//after this pass through the loop, always do the else clause
	saved_edge=next_edge-1;
      }else{
	//connect right side of previously created edge
	edges[saved_edge].ref=h;
	edges[saved_edge].ret=next_edge;
      }
      edges[next_edge].to=i;
      edges[next_edge].from=edges[h].from;
      edges[h].ret=saved_edge;
      edges[h].ref=next_edge;
      edges[next_edge].lef=h;
      edges[next_edge].let=saved_edge;
      //swap it if it needs it
      edges[h].swap(h,edges,p);
      //keep the edge for below - this is necessary in case the upper hull wasn't visible
      saved_edge=next_edge;
      next_edge++;
      //delete lower edge from the hull
      //and go round the hull in the negative direction
      lower_hull_pos=hull.delNextNeg(lower_hull_pos);
    }
    //add edge to the hull
    //set the new lower edge to the most recently created lower edge
    //if no lower hull pos was visible, use the saved edge from the upper pos earlier on
    lower_hull_pos=hull.addAfter(lower_hull_pos,saved_edge);
  }

  { 
    int i=0;
    while(edges[i].from != -1 ) i++;
    *pnedges = i;
  }
  *edges_ret = edges;

}
 
#include "heapsort.h"

void sort_triangulate(int npoints, point *p, int *pnedges, edge** edges_ret){

  {
    for(int ip=0;ip!=npoints;++ip)
      p[ip].id = ip;
  }

  //sort the points - note that the point class defines the
  // < operator so that the sort is on the x co-ordinate
  //array p will be replaced with the array sorted in x
  heapsort(npoints,p);

#if defined(TIMING)
  {
    time_t t1 = time(NULL);
    clock_t tick1 = clock();
#endif

    //triangulate the set of points
    triangulate(npoints,p,pnedges, edges_ret);


#if defined(TIMING)
    time_t t2 = time(NULL);
    clock_t tick2 = clock();
    cout << "elapsed time (time) = " << difftime(t2,t1) << " s" << endl;
    cout << "elapsed time (clock)= " << (double)(tick2-tick1)/CLOCKS_PER_SEC << " s" << endl;
  }
#endif
}

// Auxilary class used when building the element to node connectivity
// table
class edge_auxi_t {
public:
  edge_auxi_t()
    :
    left_visited_(false), right_visited_(false),
    ie_left(-2), ie_right(-2)
  {}
  bool left_visited() const { return left_visited_; }
  bool right_visited() const { return right_visited_; }
  int ielem_left() const;
  int ielem_right() const;
  void mark_left(int ielem);
  void mark_right(int ielem);
private:
  bool left_visited_, right_visited_;
  // element on the right and left side
  int ie_left, ie_right;
};

void edge_auxi_t::mark_left(int ielem) { 
  left_visited_ = true;
  ie_left = ielem;
}
void edge_auxi_t::mark_right(int ielem) { 
  right_visited_ = true;
  ie_right = ielem;
} 
int edge_auxi_t::ielem_left() const { assert(left_visited_); return ie_left; }
int edge_auxi_t::ielem_right() const { assert(right_visited_); return ie_right; }


// auxiliary class to get clockwise and counter clockwise edges around a node
class oriented_edge {
  int _edge;
  bool _orientation;
  int o() const { return _orientation; }
public:
  oriented_edge():
    _edge(-1),
    _orientation(true) {}
  oriented_edge(int e, bool o):
    _edge(e),
    _orientation(o) {}
  oriented_edge(const oriented_edge & _e):
    _edge(_e.e()),
    _orientation(_e.o()) {}
  int e() const { return _edge; }
  oriented_edge next_ccw_around_to(const edge* edges) const;
  oriented_edge next_cw_around_to(const edge* edges) const ;
  oriented_edge ccw_edge_around_to(const edge* edges) const;
};

oriented_edge oriented_edge::next_ccw_around_to(const edge* edges) const {
  int ires;
  bool bres = true;
  if (o()){
    ires = edges[e()].ret;
  } else {
    ires = edges[e()].lef;
  }
  if (ires != -1){
    bres = edges[ires].let == e();
    if (!bres)
      assert( edges[ires].ref == e());
  }
  return oriented_edge(ires,bres);
}

oriented_edge oriented_edge::next_cw_around_to(const edge* edges) const {
  int ires;
  bool bres = true;
  if (o()){
    ires = edges[e()].let;
  } else {
    ires = edges[e()].ref;
  }
  if (ires != -1){
    bres = edges[ires].ret == e();
    if (!bres)
      assert( edges[ires].lef == e());
  }
  return oriented_edge(ires,bres);
}

// give counter clockwise edge
oriented_edge oriented_edge::ccw_edge_around_to(const edge* edges) const {
  oriented_edge ccw = next_ccw_around_to(edges);
  if (ccw.e() == -1) {
    // iterated on clockwise edges
    oriented_edge e1(*this);
    for(;;){
      const oriented_edge enext = e1.next_cw_around_to(edges);
      if (enext.e() == -1) 
	break;
      e1 = enext;
    }
    ccw = e1;
  }
  return ccw;
}


// mark as visited the side of iedge_markable that points to
// iedge_orig
static
void mark_as_visited(int iedge_markable, int iedge_orig,
		     const edge* edges, edge_auxi_t* edges_visit,
		     int ielem){
  if (edges[iedge_markable].lef == iedge_orig ||
      edges[iedge_markable].let == iedge_orig){
    edges_visit[iedge_markable].mark_left(ielem);
  } else {
    assert(edges[iedge_markable].ref == iedge_orig ||
	   edges[iedge_markable].ret == iedge_orig);
    edges_visit[iedge_markable].mark_right(ielem);
  }
}

// check geometry properties of the edges in the elems table
static
void sanity_check_elems(const edge* edges, const elem *elems, 
			const edge_auxi_t *edges_visit, int ielem){
  if (elems[ielem].eo2){
    assert(edges[elems[ielem].e2].from == elems[ielem].p2);
    assert(edges[elems[ielem].e2].to == elems[ielem].p1);
    assert(edges_visit[elems[ielem].e2].ielem_right() == ielem);
  } else {
    assert(edges[elems[ielem].e2].from == elems[ielem].p1);
    assert(edges[elems[ielem].e2].to == elems[ielem].p2);
    assert(edges_visit[elems[ielem].e2].ielem_left() == ielem);
  }
  if (elems[ielem].eo1){
    assert(edges[elems[ielem].e1].from == elems[ielem].p1);
    assert(edges[elems[ielem].e1].to == elems[ielem].p3);
    assert(edges_visit[elems[ielem].e1].ielem_right() == ielem);
  } else {
    assert(edges[elems[ielem].e1].from == elems[ielem].p3);
    assert(edges[elems[ielem].e1].to == elems[ielem].p1);
    assert(edges_visit[elems[ielem].e1].ielem_left() == ielem);
  }
  if (elems[ielem].eo3){
    assert(edges[elems[ielem].e3].from == elems[ielem].p3);
    assert(edges[elems[ielem].e3].to == elems[ielem].p2);
    assert(edges_visit[elems[ielem].e3].ielem_right() == ielem);
  } else {
    assert(edges[elems[ielem].e3].from == elems[ielem].p2);
    assert(edges[elems[ielem].e3].to == elems[ielem].p3);
    assert(edges_visit[elems[ielem].e3].ielem_left() == ielem);
  }
}

// build connectivity table element to node and edges
// build pelem, pelems_ret
void build_elem_table(int npoints, const point *p, int nedges, const edge* edges,
		      int *pnelem, elem** pelems_ret){
  // Euler invariant for a triangulation:
  const int nelem = 1 + nedges - npoints;
  elem *elems = new elem[nelem]; 
  edge_auxi_t *edges_visit = new edge_auxi_t[nedges];
  { 
    // build non oriented edges for each element
    // scan all edges. Assign element to right and left
    // and mark the visited sides for each edges

    // build edges per element
    int ielem = 0;
    for(int iedge=0;iedge<nedges;iedge++) {
      // left
      if (! edges_visit[iedge].left_visited()) {
	int ielem_current = -1;
	if (edges[iedge].lef == -1) {
	  assert(edges[iedge].let == -1);
	} else {
	  // don't bother with orientation at the moment
	  ielem_current = ielem;
	  elems[ielem].e1 = iedge;
	  elems[ielem].e2 = edges[iedge].lef;
	  mark_as_visited(edges[iedge].lef,iedge,edges, edges_visit, ielem_current);
	  elems[ielem].e3 = edges[iedge].let;
	  mark_as_visited(edges[iedge].let,iedge,edges, edges_visit, ielem_current);
	  ielem++;
	}
	edges_visit[iedge].mark_left(ielem_current);
      }
      // right
      if (! edges_visit[iedge].right_visited()) {
	int ielem_current = -1;
	if (edges[iedge].ref == -1) {
	  assert(edges[iedge].ref == -1);
	} else {
	  // don't bother with orientation at the moment
	  ielem_current = ielem;
	  elems[ielem].e1 = iedge;
	  elems[ielem].e2 = edges[iedge].ref;
	  mark_as_visited(edges[iedge].ref,iedge,edges, edges_visit, ielem_current);
	  elems[ielem].e3 = edges[iedge].ret;
	  mark_as_visited(edges[iedge].ret,iedge,edges, edges_visit, ielem_current);
	  ielem++;
	}
	edges_visit[iedge].mark_right(ielem_current);
      }
    }
    assert(ielem == nelem);
  }
  for(int ielem=0;ielem<nelem;ielem++){
    {
      // build vertices per element
#define SWAP_E(EI,EJ) \
          { int etemp = elems[ielem].EI; \
	  elems[ielem].EI = elems[ielem].EJ; \
	  elems[ielem].EJ = etemp; }
      
      elems[ielem].p1 = edges[elems[ielem].e2].to;
      elems[ielem].p2 = edges[elems[ielem].e2].from;
      // using e3
      int itemp = edges[elems[ielem].e3].from;
      if (itemp == elems[ielem].p1 || itemp == elems[ielem].p2){
	elems[ielem].p3 = edges[elems[ielem].e3].to;
	if (itemp == elems[ielem].p1) { // e3 joins p1-p3
          SWAP_E(e1,e3);
	}
      } else {
	elems[ielem].p3 = itemp;
	if (edges[elems[ielem].e3].to == elems[ielem].p1){ // e3 joins p1-p3
          SWAP_E(e1,e3);
	}
      }
      // orientation: counter clockwise
      // angle(p2p1.p2p3) must be negative (counter clockwise)
      // <=> sin(p2p1.p2p3) <= 0 <=> vect_prod(p2p1,p2p3) <=0
      const double v=
	vecprod(elems[ielem].p2,elems[ielem].p1,elems[ielem].p3,p);
      // points should not be aligned.
      assert(v != 0.);
      if (v>0) {
	SWAP_E(p1,p3)
	  SWAP_E(e2,e3);
      }
      // edges orientations
      elems[ielem].eo2 = (edges[elems[ielem].e2].from == elems[ielem].p2);
      elems[ielem].eo3 = (edges[elems[ielem].e3].from == elems[ielem].p3);
      elems[ielem].eo1 = (edges[elems[ielem].e1].from == elems[ielem].p1);
#undef SWAP_E
    }
#if 1
    sanity_check_elems(edges, elems, edges_visit, ielem);
#endif
    {
      // build t1, t2, t3
      if (elems[ielem].eo2){
	elems[ielem].t3 = edges_visit[elems[ielem].e2].ielem_left();
      } else {
	elems[ielem].t3 = edges_visit[elems[ielem].e2].ielem_right();
      }
      if (elems[ielem].eo1){
	elems[ielem].t2 = edges_visit[elems[ielem].e1].ielem_left();
      } else {
	elems[ielem].t2 = edges_visit[elems[ielem].e1].ielem_right();
      }
      if (elems[ielem].eo3){
	elems[ielem].t1 = edges_visit[elems[ielem].e3].ielem_left();
      } else {
	elems[ielem].t1 = edges_visit[elems[ielem].e3].ielem_right();
      }

    }
  }
#if 0
  {
    for(int ielem=0;ielem<nelem;ielem++){
      cout << "elem=" << ielem
	   << " p1=" << elems[ielem].p1 << " (" << p[elems[ielem].p1].x <<"," << p[elems[ielem].p1].y << ")"
	   << " p2=" << elems[ielem].p2 << " (" << p[elems[ielem].p2].x <<"," << p[elems[ielem].p2].y << ")"
	   << " p3=" << elems[ielem].p3 << " (" << p[elems[ielem].p3].x <<"," << p[elems[ielem].p3].y << ")"
	   << endl;
    }
  }
#endif
  {
    // test ccw edge code
    for(int iedge=0; iedge<nedges; iedge++){
      const oriented_edge e1(iedge,true);
      const oriented_edge ccw_to = e1.ccw_edge_around_to(edges);
      const oriented_edge e2(iedge,false);
      const oriented_edge ccw_from = e2.ccw_edge_around_to(edges);
#if 0
      cout << "edge=" << iedge
	   << " ret=" << edges[iedge].ret 
	   << " lef=" << edges[iedge].lef
	   << " ccw_to=" << ccw_to.e()
	   << " ccw_from=" << ccw_from.e()
	   << endl;
#endif
    }
  }

  delete [] edges_visit; edges_visit = NULL;
  *pnelem = nelem;
  *pelems_ret = elems;
}


void sort_triangulate(int npoints, point *p,
		      int *pnedges, edge** edges_ret,
		      int *pnelem, elem** pelems_ret){
  sort_triangulate(npoints,p,pnedges, edges_ret);
  build_elem_table(npoints, p, *pnedges, *edges_ret, pnelem, pelems_ret);
}
