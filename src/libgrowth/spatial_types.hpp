#ifndef SPATIAL_TYPES_H
#define SPATIAL_TYPES_H

// C++ include
#include <cmath>
#include <memory>
#include <vector>

// GEOS includes
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/prep/PreparedGeometry.h>
#include <geos/index/strtree/STRtree.h>
#include <geos/operation/buffer/BufferParameters.h>
// define the proper types to use in shared_ptr (avoid memory leak)
#define GEOSGeometry geos::geom::Geometry
#define GEOSPreparedGeometry geos::geom::prep::PreparedGeometry
#define GEOSCoordSequence geos::geom::CoordinateSequence
#define GEOSSTRtree geos::index::strtree::STRtree
#define GEOSBufferParams geos::operation::buffer::BufferParameters
#define GEOSSTRtree_t GEOSSTRtree
#define GEOSBufParams_t GEOSBufferParams
// include C-API
#include <geos_c.h>
// define missing types
typedef GEOSPreparedGeometry GEOSPrepGeom_t;
typedef GEOSGeometry GEOSGeometry_t;


namespace growth
{

/*
 * GEOS pointers
 */

typedef std::shared_ptr<GEOSGeometry_t> GeomPtr;
typedef std::shared_ptr<GEOSPrepGeom_t> PrepGeomPtr;


/*
 * Geometrical objects
 */

typedef std::vector<double> Vector;


class Point
{
  public:
    Point();
    Point(double x, double y);
    Point(const Point &pt);

    double operator[](const int idx);
    double at(const int idx) const;

  private:
    double x_;
    double y_;
};


/*
 * Growth cone motion
 */

typedef struct Filopodia
{
    std::vector<double> directions;
    std::vector<double> normal_weights;
    int size;
    double wall_affinity;
    double finger_length;
} Filopodia;


typedef struct Move
{
    Point position;
    double module;
    double speed;
    double angle;
    double sigma_angle;
    Move()
        : position(Point())
        , module(0)
        , speed(0)
        , angle(0)
        , sigma_angle(0)
    {
    }
} Move;


/*
 * Tools
 */

inline double _rad_from_deg(double deg_angle)
{
    return deg_angle / 180. * M_PI;
}


inline double _deg_from_rad(double rad_angle)
{
    return rad_angle * 180. / M_PI;
}


/*
 * Shape
 */

class Shape
{
    Shape();
};
}

#endif /* SPATIAL_TYPES_H */
