//
//  AreaCalculator.h
//  gmetric
//
//  Created by Leonhard Spiegelberg on 18.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef __gmetric__AreaCalculator__
#define __gmetric__AreaCalculator__

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

// GNU Scientific Library
#include <gsl/gsl_integration.h>

#include "NCGrid.h"
#include "Vector2.h"
#include "Vector3.h"
#include "MathFunctions.h"
#include "helper.h"

/// small enum for sides
enum CellSide {
    SIDE_TOP = 0,
    SIDE_BOTTOM = 1,
    SIDE_LEFT = 2,
    SIDE_RIGHT = 3,
    SIDE_UNKNOWN = 4,
    SIDE_TRIANGLE = 5 // for volume calculation!
};

/// struct to hold information for intersection point
template<typename T> struct SIntersectionPoint {
    Vector2<T> point;
    CellSide side;
};

/// struct to hold two intersection points(forming line segment)
/// an assigned depth value
/// and a distance(could be norm2 or calculated with gsl)
template<typename T> struct SLineSegment{
    Vector2<T>  X;
    Vector2<T>  Y;
    T           depth;
    T           distance;
};

/// class to calculate the area metric (define it as singleton)
class AreaCalculator {
private:
    /// GSL workspace
    gsl_integration_workspace * workspace;
    
    
    /// constructor allocating space for GSL integration space
    AreaCalculator() {
        // GSL, allocate integration space to solve real world integrals
        workspace= gsl_integration_workspace_alloc(3000);
    }
    
    // static function to calculate real world data
    // params should be an array of three double's (for a,b,c of a line)
    static double df (double x, void * params) {
        
        double *coeff = (double*)params;
        
        double t = 0;
        
        if(abs(coeff[1]) >= EPSILON) {
            // case \tilde{b} != 0
            assert(coeff[1] != 0.0);
            
            double a = -coeff[0] / coeff[1];
            double b = -coeff[2] / coeff[1];
            
            double cosaxb = cos(a * x + b);
            double sinaxb = sin(a * x + b);
            double u1 = -sin(x) * cosaxb - a * sinaxb * cos(x);
            double u2 =  cos(x) * cosaxb - a * sinaxb * sin(x);
            double u3 =  a * cosaxb;
            t = u1 * u1 + u2 * u2 + u3 * u3;
        }
        else {
            
            // case \tilde{b} == 0
            assert(coeff[0] != 0.0);
            
            double c = -coeff[2] / coeff[0];
            
            double u1 = -cos(c) * sin(t);
            double u2 = -sin(c) * sin(t);
            double u3 =  cos(x);
            t = u1 * u1 + u2 * u2 + u3 * u3;
        }
        
        return sqrt(t);
    }
    
    /// helper function, to find for a grid with arbitrary coordinates an index bounding box for given bounding box by min/max
    void findBoundingBoxIndices(NCGrid *grid, const double2 min, const double2 max, int& si, int & ei, int& sj, int & ej);
    
    bool    splineLineIntersectionLon(const double2 A, const double2 B, const double *spline, double splineALon, double splineBLon, double3 *out, unsigned int *count);
    bool    splineLineIntersectionLat(const double2& A, const double2& B, const double *spline, const          double splineALat, const double splineBLat, double3 *out, unsigned int *count);
    
    std::vector<SLineSegment<double> > splineSegmentsInsideCell(const double2& A, const double2& B, const double2& topleft, const double2& topright, const double2& bottomright, const double2& bottomleft, double *xTopSpline, double *xBottomSpline, double *yLeftSpline, double *yRightSpline);
    
    std::string printOutSplineCode(const double *a, double start, double end, char mode, const double3& col = double3(0.7, 0.7, 0.7), const double linewidth=1.0);
    bool llPointInSplineCell(const double2& p, const double *left, const double * right, const double *top, const double *bottom);
    double splineLineIntegral(const double a, const double b, const double *spline, const double2& A, const double2& B);
    
    /// helper function for sorting SIntersectionPoints after lon/lat
    // return values:
    // <0	The element pointed by p1 goes before the element pointed by p2
    //  0	The element pointed by p1 is equivalent to the element pointed by p2
    // >0	The element pointed by p1 goes after the element pointed by p2
    template<typename T> static int compIntersectionPoints(const void *_a, const void *_b) {
        SIntersectionPoint<T> *a = (SIntersectionPoint<T>*)_a;
        SIntersectionPoint<T> *b = (SIntersectionPoint<T>*)_b;
        if(a->point.lon < b->point.lon)return -1;
        else {
            if(a->point.lon > b->point.lon)return 1;
            else {
                if(a->point.lat < b->point.lat)return -1;
                else if(a->point.lat > b->point.lat)return 1;
                else return 0;
            }
        }
    }
    
    /// helper function for sorting double3 after lon/lat
    /// important: double3 should have z = 1!
    // return values:
    // <0	The element pointed by p1 goes before the element pointed by p2
    //  0	The element pointed by p1 is equivalent to the element pointed by p2
    // >0	The element pointed by p1 goes after the element pointed by p2
    static int compdouble3Points(const void *_a, const void *_b) {
        double3 a = *((double3*)_a);
        double3 b = *((double3*)_b);
        
        assert(fequals(a.z, 1));
        assert(fequals(b.z, 1));
        
        if(a.x < b.x)return -1;
        else {
            if(a.x > b.x)return 1;
            else {
                if(a.y < b.y)return -1;
                else if(a.y > b.y)return 1;
                else return 0;
            }
        }
    }
    
    
    /// helper function to check if line points occur in any cell
    bool pointsInAnyLinearCell(NCGrid *grid, const int si, const int ei, const int sj, const int ej, const double2& A, const double2& B);
    
    /// helper function to check if line points occur in any cell
    bool pointsInAnySplineCell(NCGrid *grid, const int si, const int ei, const int sj, const int ej, const double2& A, const double2& B);
    
public:
    ~AreaCalculator() {
        // free GSL space
        gsl_integration_workspace_free(workspace);
    }
    
    /// singleton getter
    static AreaCalculator& instance() {
        static AreaCalculator ac;
        return ac;
    }
    
    /// calculate area, given by a NCGrid, where edges are set to lines
    /// empty string as filename means no output
    double  linear(NCGrid *grid, const double2 _A, const double2 _B, const std::string filename, const bool realWorldData, const double R = 1.0);
    
    /// calculate area, given by a NCGrid, with edges as splines
    /// empty string as filename means no output
    double  spline(NCGrid *grid, const double2 _A, const double2 _B, const std::string filename, const bool realWorldData, const double R = 1.0);
    
    /// helper function to calculate distance for two longitude/latitude points either on a sphere or just per Pythagorean norm
    double llDistanceUnitSphere(const double2 A, const double2 B, bool realWorldData);
    
    friend class AreaCalculatorTest;
    
    // avoid writing functions doubled, volume calculator is a friend!
    friend class VolumeCalculator;
};

#endif /* defined(__gmetric__AreaCalculator__) */
