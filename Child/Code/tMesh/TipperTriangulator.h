//-*-c++-*-

/***************************************************************************/
/**
**  @file
**  @brief definitions for Tipper triangulation
*/
/***************************************************************************/

#include <iosfwd>

#ifndef BOOL
#define BOOL(x) (x)
#endif
#if defined(__SUNPRO_CC)
# if __SUNPRO_CC==0x420 && !defined(ENUM_BOOL_DEFINED)
#  define ENUM_BOOL_DEFINED 1
typedef enum { false=0, true } bool;
#undef BOOL
#define BOOL(x) ((x)?true:false)
# endif
#endif

#define DEBUG_PRINT 1

class point;
class edge;
class elem;
class oriented_edge;
void tt_sort_triangulate(int npoints, point *p,
			 int *pnpoints_unique,
			 int *pnedges, edge** edges_ret);
void tt_sort_triangulate(int npoints, point *p,
			 int *pnpoints_unique,
			 int *pnedges, edge** edges_ret,
			 int *pnelem, elem** pelems_ret);
void tt_build_elem_table(int npoints, const point *p,
			 int nedges, const edge* edges,
			 int *pnelem, elem** pelems_ret);
void tt_build_spoke(int npoints, int nedges, const edge* edges,
		    oriented_edge** poedge);
void tt_error_handler(void);

class point{
public:
  point() : _id(-1) { _XY[0] = 0.; _XY[1] = 0.;}
  point(double ix,double iy) : _id(-1) { _XY[0] = ix; _XY[1] = iy;}
  point(double ix,double iy,int iid) : _id(iid) { _XY[0] = ix; _XY[1] = iy; }
  point(const point& p) : _id(p.id()) { _XY[0] = p.x(); _XY[1] = p.y(); }
  const point &operator=( const point &p );
  bool operator < (const point& p) const;
  double dot(const point& p) const {return (x()*p.x()+y()*p.y());}
#if defined(DEBUG_PRINT)
  void print () const;
#endif
  void write(std::ofstream& f) const;
  double x() const { return _XY[0]; }
  double y() const { return _XY[1]; }
  const double *XY() const { return _XY; }
  int id() const { return _id; }
private:
  double _XY[2]; // so that it can be accessed as an array
  int _id;
};

inline point operator-(const point& lhs, const point& rhs) {
  return point(lhs.x()-rhs.x(), lhs.y()-rhs.y());
}


class edge{
  const edge &operator=( const edge & );
  edge(const edge&);
public:
  edge():
    from(end),to(end),
    lef(none),let(none),ref(none),ret(none) {}
#if defined(DEBUG_PRINT)
  void print(const point p[]) const;
#endif
  void write(std::ofstream& f,const point p[]) const;
  bool visible(const point p[],int i) const;

  enum { none = -1, end = -2 }; // must be strictly negative
  int from,to;
  int lef,let,ref,ret;
  //                       .
  //        to             .
  //    let /|\ ret        .
  //       / | \           .
  //       \ | /           .
  //    lef \|/ ref        .
  //       from            .
};

// compile-time assertions
void require_edge_c1(int d[( edge::none<0 ?1:-1)]);
void require_edge_c2(int d[( edge::end<0 ?1:-1)]);

// auxiliary class to get clockwise and counter clockwise edges around a node
class oriented_edge {
  int _edge;
  bool _orientation;
public:
  oriented_edge():
    _edge(edge::none),
    _orientation(true) {}
  oriented_edge(int _e, bool _o):
    _edge(_e),
    _orientation(_o) {}
  oriented_edge(const oriented_edge & _e):
    _edge(_e.e()),
    _orientation(_e.o()) {}
  const oriented_edge &operator=( const oriented_edge &);
  int e() const { return _edge; }
  bool o() const { return _orientation; }
  bool nonvalid() const { return (e() == edge::none ? true:false); }
  void set(int e1, bool o1) { _edge = e1; _orientation = o1; }
  oriented_edge next_ccw_around_from(const edge* edges) const;
  oriented_edge next_cw_around_from(const edge* edges) const ;
  oriented_edge ccw_edge_around_from(const edge* edges) const;
};

// connectivity table element to node and edge
class elem {
  const elem &operator=( const elem & );
  elem(const elem&);
public:
  elem() :
    p1(-2), p2(-2), p3(-2),
    e1(-2), e2(-2), e3(-2),
    t1(-2), t2(-2), t3(-2),
    eo1(false), eo2(false), eo3(false)
  {}
  int p1, p2, p3;  // nodes
  int e1, e2, e3;  // edges
  int t1, t2, t3;  // triangles (or elements)
  bool eo1, eo2, eo3; // orientation of edges
  enum{ none = -1 }; // no neighbour triangle, must be strictly negative
  //
  //         P1       .
  //        -/\       .
  //   T3 e2/  \e1 T2 .
  //       /    \     .
  //      /      \-   .
  //     P2------P3   .
  //       | e3       .
  //         T1       .

};

// compile-time assertions
void require_elem_c1(int d[( elem::none<0 ?1:-1)]);


