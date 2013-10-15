//
//  VolumeCalculator.h
//  gmetric
//
//  Created by Leonhard Spiegelberg on 30.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef __gmetric__VolumeCalculator__
#define __gmetric__VolumeCalculator__

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

// area Calculator for helper functions
#include "AreaCalculator.h"


// helper define to plot out matlab file
#define MATLAB_FILE "/Users/leonhardspiegelberg/Documents/MATLAB/vcgrid.m"
#define MATLAB // uncomment to have a MATLAB file produced!

// used for sorting
struct SAreaPoint {
    SIntersectionPoint<double> p;
    double angle;
};

/// class to calculate the volume metric (define it as singleton)
class VolumeCalculator {
private:
    /// GSL workspace
    gsl_integration_workspace * workspace;
    
    /// debug file for matlab
#ifdef MATLAB
    std::ofstream ofs;
#endif
    
    /// constructor allocating space for GSL integration space
    VolumeCalculator() {
        // GSL, allocate integration space to solve real world integrals
        workspace= gsl_integration_workspace_alloc(8000);
        
#ifdef MATLAB
        ofs.open(MATLAB_FILE);
#endif
    }
    
#ifdef MATLAB
    void mLine(const double2& A, const double2& B, const double3& col) {
        ofs<<"plot(["<<A.x<<" "<<B.x<<"],["<<A.y<<" "<<B.y<<"], 'Color',["<<col.x<<" "<<col.y<<" "<<col.z<<"]); hold on;"<<std::endl;
    }
    void mPoint(const double2& P, double size, const double3& col) {
        ofs<<std::endl<<"scatter("<<P.x<<","<<P.y<<", "<<size<<", 'o', 'fill', 'MarkerFaceColor',["<<col.x<<" "<<col.y<<" "<<col.z<<"], 'MarkerEdgeColor', 'none'); hold on;\n";
    }
#endif
    
    /// helper function for sorting SAreaPoints after angle & side
    // return values:
    // <0	The element pointed by p1 goes before the element pointed by p2
    //  0	The element pointed by p1 is equivalent to the element pointed by p2
    // >0	The element pointed by p1 goes after the element pointed by p2
    static int compAreaPoints(const void *_a, const void *_b) {
        SAreaPoint *a = (SAreaPoint*)_a;
        SAreaPoint *b = (SAreaPoint*)_b;
        // get order
        int orda = sideOrder(a->p.side);
        int ordb = sideOrder(b->p.side);
        if(orda < ordb) {
            return -1;
        }
        else if(orda > ordb) {
            return 1;
        }
        else {
            // sort according to side after lon/lat!
            // notice that here sides are equal!
            switch(a->p.side) {
                case SIDE_TOP:
                    if(a->p.point.lon < b->p.point.lon)return -1;
                    else if(a->p.point.lon > b->p.point.lon)return 1;
                    else return 0;
                    break;
                case SIDE_RIGHT:
                    if(a->p.point.lat > b->p.point.lat)return -1;
                    else if(a->p.point.lat < b->p.point.lat)return 1;
                    else return 0;
                    break;
                case SIDE_BOTTOM:
                    if(a->p.point.lon > b->p.point.lon)return -1;
                    else if(a->p.point.lon < b->p.point.lon)return 1;
                    else return 0;
                    break;
                case SIDE_LEFT:
                    if(a->p.point.lat < b->p.point.lat)return -1;
                    else if(a->p.point.lat > b->p.point.lat)return 1;
                    else return 0;
                    break;
                case SIDE_TRIANGLE:
                    // use here angle!
                    if(a->angle < b->angle)return -1;
                    else if(a->angle > b->angle)return 1;
                    else return 0;
                    break;
                default:
                    return 0;
                    break;
            }
            return 0;
        }
    }
    
    /// little helper for the side order
    static int sideOrder(CellSide side) {
        // order is top - right - bottom - left
        switch(side) {
            case SIDE_TOP:
                return 1;
                break;
            case SIDE_RIGHT:
                return 2;
                break;
            case SIDE_BOTTOM:
                return 3;
                break;
            case SIDE_LEFT:
                return 4;
                break;
            default:
                return 0xFFFF;
                break;
        }
    }
    
    /// static function to evaluate splineline integral over sphere's surface
    /// formula is int_l1^l2 sin(b_3 l^3 + b_2 l^2 + b_1 l + b_0) - sin(a_1 l + a_0) dl
    static double slf(double x, void *params) {
        assert(params);
        
        double *c = (double*)params;
        double b_0 = c[0];
        double b_1 = c[1];
        double b_2 = c[2];
        double b_3 = c[3];
        
        double a_0 = c[4];
        double a_1 = c[5];
        
        double xx = x*x;
        double xxx = x * xx;
        
        return sin(b_0 + b_1 * x + b_2 * xx + b_3 * xxx) - sin(a_0 + a_1 * x);
    }
    
    
    /// helper function to calculate area of a (convex) polygon
    /// center is the reference point from which angles shall be determined
    double  llPolygonAreaUnitSphere(double2 *points, const int count, const double2& center, const bool realWorldData);
    
    /// helper function to calculate area of a triangle
    double llTriangleUnitSphere(const double2& _A, const double2& _B, const double2& _C, const bool realWorldData);
    
    
    /// helper to calc area of a triangle with a spline segment
    /// the given spline is function of lambda! splineA & splineB are the two points of the spline segment and C another point!
    double llTriangleWithSplineLon(const double2& splineA, const double2& splineB, const double2& C, const double *spline, const bool realWorldData);
    
    /// helper to calc area of a triangle with a spline segment
    /// the given spline is function of phi! splineA & splineB are the two points of the spline segment and C another point!
    double llTriangleWithSplineLat(const double2& splineA, const double2& splineB, const double2& C, const double *spline, const bool realWorldData);
    
    
    /// helper function to intersect all spline sides with all cell sides
    void intersectTriangleSidesWithCellSides(const double2& A, const double2& B, const double2& C, const double2 *p, const double * topSpline, const double *bottomSpline, const double *leftSpline, const double *rightSpline, SIntersectionPoint<double> *out, int *count);
    
    /// helper function to intersect line with all splines
    void intersectLineSegmentWithCellSides(const double2& _A, const double2& _B, const double2 *p, const double * topSpline, const double *bottomSpline, const double *leftSpline, const double *rightSpline, SIntersectionPoint<double> *out, int *count);
    
    /// helper function to get area of the polygon where borders can be splines
    /// center is the reference point from which angles shall be determined
    double  llPolygonWithSplinesAreaUnitSphere(SIntersectionPoint<double> *points, const int count, const double2& center, const double * topSpline, const double *bottomSpline, const double *leftSpline, const double *rightSpline, const bool realWorldData);
    
    /// helper function to integrate over the region between spline and line
    /// input values are in degrees!
    double splineLineIntegral(const double _a, const double _b, const double *spline, const double2& _A, const double2& _B, const bool realWorldData);

public:
    ~VolumeCalculator() {
        // free GSL space
        gsl_integration_workspace_free(workspace);
    }
    
    /// singleton getter
    static VolumeCalculator& instance() {
        static VolumeCalculator vc;
        return vc;
    }
    
    /// calculate volume of triangle, given by a NCGrid, where edges are set to lines
    /// empty string as filename means no output
    double  linear(NCGrid *grid, const double2 _A, const double2 _B, const double2 _C, const std::string filename, const bool realWorldData, const double R = 1.0);
    
    /// calculate volume, given by a NCGrid, with edges as splines
    /// empty string as filename means no output
    double  spline(NCGrid *grid, const double2 _A, const double2 _B, const double2 _C, const std::string filename, const bool realWorldData, const double R = 1.0);
    
    friend class VolumeCalculatorTest;
};


#endif /* defined(__gmetric__VolumeCalculator__) */
