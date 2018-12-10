#ifndef SPATIAL_TYPES_H
#define SPATIAL_TYPES_H
#ifndef GEOS_USE_ONLY_R_API
#define GEOS_USE_ONLY_R_API
#endif

// C++ include
#include <cmath>
#include <memory>
#include <vector>
#include <unordered_map>
#include <tuple>

// Boost includes
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/variant.hpp>

// include GEOS C-API
#include <geos_c.h>


namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;


namespace growth
{

/*
 * Boost objects
 */

typedef bg::model::d2::point_xy<double> BPoint;
typedef bg::model::multi_point<BPoint> BMultiPoint;
typedef bg::model::segment<BPoint> BSegment;
typedef bg::model::linestring<BPoint> BLineString;
typedef bg::model::multi_linestring<BLineString> BMultiLineString;
typedef bg::model::ring<BPoint> BRing;
typedef bg::model::polygon<BPoint> BPolygon;
typedef bg::model::multi_polygon<BPolygon> BMultiPolygon;
typedef bg::model::box<BPoint> BBox;

typedef boost::variant<BPolygon, BMultiPolygon> BGeometry;

typedef std::shared_ptr<BPolygon> BPolygonPtr;
typedef std::shared_ptr<BMultiPolygon> BMultiPolygonPtr;
typedef std::shared_ptr<BGeometry> BGeometryPtr;


typedef std::vector<BPolygonPtr>::const_iterator seg_it;
typedef boost::iterator_range< seg_it > seg_range;


/*
 * Smart pointers
 */

class Area;
typedef std::shared_ptr<Area> AreaPtr;


/*
 * Geometrical objects
 */

typedef std::vector<double> Vector;
typedef std::array<double, 3> PointArray;
typedef std::array<Vector, 3> PointsArray;


/*
 * Growth cone motion
 */

typedef struct Filopodia
{
    std::vector<double> directions;
    std::vector<double> normal_weights;
    int size;
    double finger_length;
    double substrate_affinity;
    double wall_affinity;
} Filopodia;


typedef struct Move
{
    BPoint position;
    double module;
    double speed;
    double angle;
    double sigma_angle;
    Move()
        : position(BPoint())
        , module(0)
        , speed(0)
        , angle(0)
        , sigma_angle(0)
    {
    }
} Move;


/*
 * R-tree and ppatial object description into tuple and hash function
 */

// ObjectInfo tuple (neuron ID, neurite name, node ID, index on branch)
typedef std::tuple<size_t, std::string, size_t, size_t> ObjectInfo;


// R-tree value
typedef std::pair<BBox, ObjectInfo> RtreeValue;


// tool to hash the ObjecInfo/Polygon map
namespace hash_tuple{

    template <typename TT>
    struct hash
    {
        size_t operator()(TT const& tt) const
        {                                              
            return std::hash<TT>()(tt);                                 
        }                                              
    };

    namespace
    {

    template <class T>
    inline void hash_combine(std::size_t& seed, T const& v)
    {
        seed ^= hash_tuple::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
            hash_combine(seed, std::get<Index>(tuple));
          }
        };

        template <class Tuple>
        struct HashValueImpl<Tuple, 0>
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            hash_combine(seed, std::get<0>(tuple));
          }
        };

    } // empty namespace

    template <typename ... TT>
    struct hash<std::tuple<TT...>> 
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {                                              
            size_t seed = 0;                             
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);    
            return seed;                                 
        }                                              
    };

} // namespace has_tuple


typedef std::unordered_map<ObjectInfo, BPolygonPtr, hash_tuple::hash<ObjectInfo>> space_tree_map;
typedef std::tuple<ObjectInfo, BBox, bool> box_tree_tuple;

} // namespace growth

#endif /* SPATIAL_TYPES_H */
