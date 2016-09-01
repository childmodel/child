//-*-c++-*- 

/***************************************************************************/
/**
**  @file geometry.h
**  @brief Definitions of some simple geometry classes.
**
**  $Id: geometry.h,v 1.6 2003-05-15 16:02:55 childcvs Exp $
*/
/***************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

class Point2D
{
  public:
   Point2D();
   Point2D( double, double );
   const Point2D &operator=( const Point2D & );
   double x;
   double y;
};

class Point3D
{
  public:
   Point3D();
   Point3D( double, double, double );
   const Point3D &operator=( const Point3D & );
   double x;
   double y;
   double z;
};

inline Point2D::Point2D() :
  x(0.), y(0.)
{}


inline Point2D::Point2D( double ix, double iy ) :
   x(ix), y(iy)
{}

inline const Point2D &Point2D::operator=(const Point2D &right )
{
  if ( &right != this ) {
    x = right.x;
    y = right.y;
  }
  return *this;
}

inline Point3D::Point3D() :
  x(0.), y(0.), z(0.)
{}

inline Point3D::Point3D( double ix, double iy, double iz ) :
  x(ix), y(iy), z(iz)
{}

inline const Point3D &Point3D::operator=(const Point3D &right )
{
  if ( &right != this ) {
    x = right.x;
    y = right.y;
    z = right.z;
  }
  return *this;
}



#endif
