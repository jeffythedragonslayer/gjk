#include "stdafx.h"
#include "Vec3.h"
#include "GJK.h"
#include "Plane.h"
using namespace std;

// Boolean Gilbert-Johnson-Keerthi Algorithm
//
// Return whether two convex objects intersect
//
// calls support() function with a search direction
// support() should be implemented as support1(dir) - support2(-dir) for the two objects
//
// Start with an arbitrary search direction, such as VEC3_UP
// Loop quits after max_iterations, most object pairs take 0-5 iterations to enclose the origin
// if something goes wrong it is safer to return true since GJK is used as a broad-phase technique on convex hulls
bool gjk_test(Vec3& direction, uint max_iterations) 
{
        // single point case is trivial
        Vec3 b = support(direction);
        if( b.is_zero() ) return true; // something went wrong, so say objects intersecting by default

        // line segment case is unrolled from the loop because once its done, we never return to it
        Vec3 a = support(-b);
        if( a.is_zero() || a.dot(b) > 0 ) return false; // something went wrong, or a is on the same side of the origin as b with respect to direction

        Vec3 ab = b - a;
        direction = ab.cross(-a).cross(ab);
        Vec3 simplex[4] = {b, a}; // initialize the simplex with 2 points
        uint dimensions = 1; // # of dimensions of the simplex (# of points == dimensions+1)
        uint i = max_iterations;

        do {
                if( direction.is_zero() ) return true; // something went wrong, say objects intersecting by default
                a = support(direction); // find new point of in minkowski difference farthest along search direction

                if( a.dot(direction) < 0 ) return false; // farthest supporting point is not past origin with respect to search direction, so origin cannot be enclosed
                simplex[++dimensions] = a; // add new point to simplex

                // at this point dimensions is either 2 or 3 (simplex is either a triangle or tetrahedron)
                // determine which voronoi region the origin is in, and update the simplex to be the feature closest to it along with search direction

                if( 2 == dimensions ){ // triangle case
                        a        = -simplex[2];
                        Vec3 ac  =  simplex[0] + a;
                        Vec3 ab  =  simplex[1] + a;
                        Vec3 abc =  ab.cross(ac);
                        if( abc.cross(ac).dot(a) > 0 ){
                                if( ac.dot(a) > 0 ){ // origin is in region AC
                                        simplex[1] = simp[2]; 
                                        direction  = ac.cross(a).cross(ac);
                                        dimensions = 1;
                                } else { // origin is in region AB
                                        simplex[0] = simp[1]; 
                                        simplex[1] = simp[2];
                                        direction  = ab.cross(a).cross(ab);
                                        dimensions = 1;
                                }
                        } else if( ab.cross(abc).dot(a) > 0 ) {
                                simplex[0] = simp[1]; // origin is in region AB
                                simplex[1] = simp[2];
                                direction = ab.cross(a).cross(ab);
                                dimensions = 1;
                        }
                        else if( abc.dot(a) > 0 ){
                                direction = abc; // origin is in region ABC, above triangle
                        } else {
                                swap(simplex[1], simp[2]); // origin is in region ABC, below triangle
                                direction = -abc;
                        }
                } else { // tetrahedron case, dimensions == 3
                        Vec3 ab = simplex[2] - simplex[3];
                        Vec3 ac = simplex[1] - simplex[3];
                        Vec3 ad = simplex[0] - simplex[3];
                        Vec3 abc = ab.cross(ac);
                        Vec3 acd = ac.cross(ad);
                        Vec3 adb = ad.cross(ab);
                        bool above_acd = point_above_plane(VEC3_ZERO, simplex[3], simplex[1], simplex[0]);
                        bool above_adb = point_above_plane(VEC3_ZERO, simplex[3], simplex[0], simplex[2]);
                        if( point_above_plane(VEC3_ZERO, simplex[3], simplex[2], simplex[1]) ){
                                if( above_acd ){
                                        simplex[0] = simp[1]; // origin is in region AC
                                        simplex[1] = simp[3];
                                        direction = ac;
                                        dimensions = 1;
                                } else if( above_adb ){
                                        simplex[0] = simp[2]; // origin is in region AB
                                        simplex[1] = simp[3];
                                        direction = ab;
                                        dimensions = 1;
                                } else {
                                        simplex[0] = simp[1]; // origin is in region ABC
                                        simplex[1] = simp[2];
                                        simplex[2] = simp[3];
                                        direction = abc;
                                        dimensions = 2;
                                }
                        } else if( above_acd ){
                                if( above_adb ){
                                        simplex[1] = simp[3]; // origin is in region AD
                                        direction = ad;
                                        dimensions = 1;
                                } else {
                                        simplex[2] = simp[3]; // origin is in region ACD
                                        direction = acd;
                                        dimensions = 2;
                                }
                        } else if( above_adb ){
                                simplex[1] = simp[2]; // origin is in region ADB
                                simplex[2] = simp[3];
                                direction = adb;
                                dimensions = 2;
                        } else return true; // origin is in region ABCD and has been enclosed
                }
        } while (--i);
        return true; // hit max iterations, so say objects intersect by default
}
