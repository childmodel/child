
#include <fstream.h>

class point;
class edge;
class elem;
void sort_triangulate(int npoints, point *p, int *pnedges, edge** edges_ret);
void sort_triangulate(int npoints, point *p,
		      int *pnedges, edge** edges_ret,
		      int *pnelem, elem** pelems_ret);
void build_elem_table(int npoints, const point *p, int nedges, const edge* edges,
		      int *pnelem, elem** pelems_ret);

class point{
public:
  point() : x(0.), y(0.) {}
  point(double ix,double iy) : x(ix), y(iy) {}
  point(const point& p) : x(p.x), y(p.y) {}
  const point &operator=( const point &p );
  int operator < (const point& p) const {return x<p.x;}
  point operator - (const point& p) const {return point(x-p.x,y-p.y);}
  point operator + (const point& p) const {return point(x+p.x,y+p.y);}
  point operator / (double f) const {return point(x/f,y/f);}
  double dot(const point& p) const {return (x*p.x+y*p.y);}
#if defined(DEBUG_PRINT)
  void print () const;
#endif
  void write(ofstream& f) const;
public:
  double x,y;
};

class edge{
  const edge &operator=( const edge & );  // assignment operator
public:
  edge(): from(-1),to(-1),lef(-1),let(-1),ref(-1),ret(-1) {}
#if defined(DEBUG_PRINT)
  void print(const point p[]) const;
#endif
  void write(ofstream& f,const point p[]) const;
  bool visible(const point p[],int i) const;
  int swap(int tint,edge e[],const point p[]);
public:
  int from,to;
  int lef,let,ref,ret;
};

class elem {
public:
  int p1, p2, p3;
  int e1, e2, e3;
};

