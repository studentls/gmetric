//
//  Metrics.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 01.08.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef __nemoExtractor__Metrics__
#define __nemoExtractor__Metrics__

#include <fstream>

#include "helper.h"
#include "Vector2.h"
#include "Vector3.h"
#include "TArray.h"
#include "MathFunctions.h"
#include "NCGrid.h"
#include "AreaCalculator.h"
#define BUFFER_SIZE 10 // for security 10 places(6 only needed)!


///// small enum for sides
//enum CellSide {
//    SIDE_TOP = 0,
//    SIDE_BOTTOM = 1,
//    SIDE_LEFT = 2,
//    SIDE_RIGHT = 3,
//    SIDE_UNKNOWN = 4
//};
//
///// struct to hold information for intersection point
//template<typename T> struct SIntersectionPoint {
//    Vector2<T> point;
//    CellSide side;
//};
//
///// struct to hold two intersection points(forming line segment)
//template<typename T> struct SLineSegment{
//    Vector2<T> X;
//    Vector2<T> Y;
//};

///// helper function for sorting SIntersectionPoints after lon/lat
//// return values:
//// <0	The element pointed by p1 goes before the element pointed by p2
////  0	The element pointed by p1 is equivalent to the element pointed by p2
//// >0	The element pointed by p1 goes after the element pointed by p2
//template<typename T> int compIntersectionPoints(const void *_a, const void *_b) {
//    SIntersectionPoint<T> *a = (SIntersectionPoint<T>*)_a;
//    SIntersectionPoint<T> *b = (SIntersectionPoint<T>*)_b;
//    if(a->point.lon < b->point.lon)return -1;
//    else {
//        if(a->point.lon > b->point.lon)return 1;
//        else {
//            if(a->point.lat < b->point.lat)return -1;
//            else if(a->point.lat > b->point.lat)return 1;
//            else return 0;
//        }
//    }
//}

class Metrics {
public:
    
    /// helper to check orientation of array
    template<typename T> static int orientationOfCoords(Vector2<T> *points, const int w, const int h) {
        
        using namespace std;
        // check if x-axis is ascending/descending
        bool errX = false;
        bool xLonAsc = false;
        bool xLonDesc = false;
        bool xLatAsc = false;
        bool xLatDesc = false;
        for(int j = 0; j < h; j++) { // go through all y layers
            for(int i = 0; i < w - 1; i++) {
                Vector2<T> p1 = points[i + j * w];
                Vector2<T> p2 = points[i +1+ j * w];
                if(points[i + j * w].lon < points[(i + 1) + j * w].lon)xLonAsc = true;
                if(points[i + j * w].lon > points[(i + 1) + j * w].lon)xLonDesc = true;
                if(points[i + j * w].lat < points[(i + 1) + j * w].lat)xLatAsc = true;
                if(points[i + j * w].lat > points[(i + 1) + j * w].lat)xLatDesc = true;
                
                
            }
        }
        
        //errX = xAsc == xDesc || ;
        
        // check if y-axis is ascending/descending
        bool errY = false;
        bool yLonAsc = false;
        bool yLonDesc = false;
        bool yLatAsc = false;
        bool yLatDesc = false;
        for(int j = 0; j < h - 1; j++) {
            for(int i = 0; i < w; i++) { // go through all x layers
                if(points[i + j * w].lon < points[i + (j + 1) * w].lon)yLonAsc = true;
                if(points[i + j * w].lon > points[i + (j + 1) * w].lon)yLonDesc = true;
                
                Vector2<T> p1 = points[i + j * w];
                Vector2<T> p2 = points[i + (j+1) * w];
                if(points[i + j * w].lat < points[i + (j + 1) * w].lat)yLatAsc = true;
                if(points[i + j * w].lat > points[i + (j + 1) * w].lat)yLatDesc = true;
            }
        }
        
        // no orientation, if asc&desc values are both true for any lon/lat for x or y
        if((xLonAsc && xLonDesc) || (yLatAsc && yLatDesc))return OM_NONE;
        
        //        // for strict orientation!
        //        // only two boolean expression are allowed to be set
        //        unsigned char sum = (unsigned char)xLonAsc + (unsigned char)xLonDesc + (unsigned char)xLatAsc + (unsigned char)xLatDesc +
        //        (unsigned char)yLonAsc + (unsigned char)yLonDesc + (unsigned char)yLatAsc + (unsigned char)yLatDesc;
        //        if(sum != 2)return OM_NONE;
        
        
        //errY = yAsc == yDesc && yAsc;
//        // primitives, created by rot90°
//        if(xLonAsc && yLatAsc)return OM_NEMO;
//        if(xLatAsc && yLonDesc)return OM_NEMO90CW;
//        if(xLonDesc && yLatDesc)return OM_NEMO180CW;
//        if(xLatDesc && yLonAsc)return OM_NEMO270CW;
//        
//        // primitives, after applying mirroring on x-axis
//        if(xLonDesc && yLatAsc)return OM_NEMOMX;
//        if(yLonDesc && xLatDesc)return OM_NEMO90CWMX;
//        if(xLonAsc && yLatDesc)return OM_NEMO180CWMX;
//        if(xLatAsc && yLonAsc)return OM_NEMO270CWMX;
        
        
        // primitives, created by rot90°
        if(xLonAsc && yLatAsc)return OM_NEMO;
        if(xLatAsc && yLonDesc)return OM_NEMO270CW;
        if(xLonDesc && yLatDesc)return OM_NEMO180CW;
        if(xLatDesc && yLonAsc)return OM_NEMO90CW;
        
        // primitives, after applying mirroring on x-axis
        if(xLonDesc && yLatAsc)return OM_NEMOMX;
        if(yLonDesc && xLatDesc)return OM_NEMO270CWMX;
        if(xLonAsc && yLatDesc)return OM_NEMO180CWMX;
        if(xLatAsc && yLonAsc)return OM_NEMO90CWMX;
        
        return -1;
    }
    
    static Vector2<int> NEMOIndexToOrientation(const int i, const int j, const int w, const int h, const int orientation) {
        switch(orientation) {
            case OM_NEMO:
                // just return as is
                return Vector2<int>(i, j);
                break;
            case OM_NEMOMX:
                // mirror on x axis
                return Vector2<int>(w - 1 - i, j);
                break;
            case OM_NEMO90CW:
                // return coordinates clockwise
                // this means
                // (1) tranpose
                // (2) reverse rows
                return Vector2<int>(h - j - 1, i);
                break;
            case OM_NEMO180CW:
                
                break;
        }
    }
    
    /// important helper to reorient an array, such that i, j directions correspond to lon/lat directions
    /// note that this can change w & h!
    template<typename S> static bool reorientArrayToNEMO(const int orientation, S *_points, int& w, int& h, bool output = false) {
        using namespace std;
        
        // output for user, if orientation is different from NEMO's default orientation
        switch (orientation) {
            case OM_NEMO:
                if(output)cout<<">> data properly aligned, no reorientation needed"<<endl;
                break;
            case OM_NEMO90CW:
                if(output) {
                    cout<<">> data rotated by 90°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a00    a01|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a10    a11|     |a00    a10|"<<endl;
                    cout<<"   +----------+     x----------+"<<endl;
                }
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO180CW:
                if(output) {
                    cout<<">> data rotated by 180°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a10    a00|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a11    a01|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------x"<<endl;
                }
                // could be done also faster... but anyway
                arrayRotate90CCW(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO270CW:
                if(output){
                    cout<<">> data rotated by 270°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------x"<<endl;
                    cout<<"   |a11    a10|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a01    a00|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------+"<<endl;
                }
                arrayRotate90CW(_points, w, h);
                break;
            case OM_NEMOMX:
                if(output){
                    cout<<">> data mirrored, reorientating..."<<endl;
                    cout<<"   x----------+     +----------x"<<endl;
                    cout<<"   |a11    a01|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a10    a00|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------+"<<endl;
                }
                arrayMirrorX(_points, w, h);
                break;
            case OM_NEMO90CWMX:
                if(output) {
                    cout<<">> data mirrored & rotated by 90°, reorientating..."<<endl;
                    cout<<"   x----------+     x----------+"<<endl;
                    cout<<"   |a01    a00|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a11    a10|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------+"<<endl;
                }
                arrayMirrorX(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO180CWMX:
                if(output) {
                    cout<<">> data mirrored & rotated by 180°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a00    a10|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a01    a11|     |a00    a10|"<<endl;
                    cout<<"   +----------+     x----------+"<<endl;
                }
                // could be done also faster... but anyway
                arrayMirrorX(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO270CWMX:
                if(output) {
                    cout<<">> data mirrored & rotated by 270°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a10    a11|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a00    a01|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------x"<<endl;
                }
                arrayMirrorX(_points, w, h);
                arrayRotate90CW(_points, w, h);
                break;
                
            default:
                cout<<">> error: no valid orientation given, stop execution"<<endl;
                return false;
                break;
        }
        return true;
    }
    
    
    /// important helper to reorient an array, such that i, j directions correspond to lon/lat directions
    /// note that this can change w & h!
    template<typename T> static bool reorientArray(Vector2<T> *_points, int& w, int& h, bool output = false) {
        using namespace std;
        
        // get orientation of data
        int orientation = Metrics::orientationOfCoords(_points, w, h);
        
        // output for user, if orientation is different from NEMO's default orientation
        switch (orientation) {
            case OM_NEMO:
                if(output)cout<<">> data properly aligned, no reorientation needed"<<endl;
                break;
            case OM_NEMO90CW:
                if(output) {
                    cout<<">> data rotated by 90°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a00    a01|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a10    a11|     |a00    a10|"<<endl;
                    cout<<"   +----------+     x----------+"<<endl;
                }
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO180CW:
                if(output) {
                    cout<<">> data rotated by 180°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a10    a00|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a11    a01|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------x"<<endl;
                }
                // could be done also faster... but anyway
                arrayRotate90CCW(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO270CW:
                if(output){
                    cout<<">> data rotated by 270°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------x"<<endl;
                    cout<<"   |a11    a10|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a01    a00|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------+"<<endl;
                }
                arrayRotate90CW(_points, w, h);
                break;
            case OM_NEMOMX:
                if(output){
                    cout<<">> data mirrored, reorientating..."<<endl;
                    cout<<"   x----------+     +----------x"<<endl;
                    cout<<"   |a11    a01|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a10    a00|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------+"<<endl;
                }
                arrayMirrorX(_points, w, h);
                break;
            case OM_NEMO90CWMX:
                if(output) {
                    cout<<">> data mirrored & rotated by 90°, reorientating..."<<endl;
                    cout<<"   x----------+     x----------+"<<endl;
                    cout<<"   |a01    a00|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a11    a10|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------+"<<endl;
                }
                arrayMirrorX(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO180CWMX:
                if(output) {
                    cout<<">> data mirrored & rotated by 180°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a00    a10|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a01    a11|     |a00    a10|"<<endl;
                    cout<<"   +----------+     x----------+"<<endl;
                }
                // could be done also faster... but anyway
                arrayMirrorX(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                arrayRotate90CCW(_points, w, h);
                break;
            case OM_NEMO270CWMX:
                if(output) {
                    cout<<">> data mirrored & rotated by 270°, reorientating..."<<endl;
                    cout<<"   x----------+     +----------+"<<endl;
                    cout<<"   |a10    a11|     |a01    a11|"<<endl;
                    cout<<"   |          | ->  |          |"<<endl;
                    cout<<"   |a00    a01|     |a00    a10|"<<endl;
                    cout<<"   +----------+     +----------x"<<endl;
                }
                arrayMirrorX(_points, w, h);
                arrayRotate90CW(_points, w, h);
                break;
                
            default:
                cout<<">> error: data is not oriented at all, stop execution"<<endl;
                return false;
                break;
        }
        return true;
    }
    
    /// helper, to calculate distance
    template<typename T> static T llDistance(const Vector2<T> A, const Vector2<T> B) {
        
        // first, keep it simple!
        return vector2Length(A - B);
    }
    
    /// helper to test, if point with lat, lon lies within cell bounded by splines top,left, right, bottom
    /// point on boundary, returns true!!!
    template<typename T> static bool llPointInSplineCell(const Vector2<T> p, const T *left, const T * right, const T *top, const T *bottom) {
        // note: left & right are splines with lat as parameter
        //       top & bottom with lon as parameter
        
        T leftVal = evalPolynomial(p.lat, left, 4);
        // check if point lies within border (use <= if boundary shall be excluded from test!)
        if(p.lon < leftVal)return false;
        T rightVal = evalPolynomial(p.lat, right, 4);
        if(p.lon > rightVal)return false;
        T topVal = evalPolynomial(p.lon, top, 4);
        if(p.lat > topVal)return false;
        T bottomVal = evalPolynomial(p.lon, bottom, 4);
        if(p.lat < bottomVal)return false;
        
        return true;
    }
    
    /// helper function to calculate are of a triangle with lines in lon/lat system
    template<typename T> static T llTriangle(const Vector2<T> A, const Vector2<T> B, const Vector2<T> C, const double R = 1.0) {
        // check first if points are sorted right
        assert(A.lon <= B.lon);
        assert(B.lon <= C.lon);
        
        
        // check if all variables are given in correct spaces
        assertlonlat(A);
        assertlonlat(B);
        assertlonlat(C);
        
#warning debug test!
        // debug reasons
        return triangleXYArea(A, B, C);
        
        // regard special cases
        // all longitudes are equal! => line has no area!
        if(A.lon == B.lon && A.lon == C.lon)return 0;
        //                                                   |\
        // first two lon are equal, meaning triangle of form |_\ e.g
        
        else if(A.lon == B.lon) { // maybe later use here fequals?
            
            // depending if a_11 = 0 or b_1 = 0, form integrals
            T intA = 0;
            T intB = 0;
            
            T a_11 = (C.lat - B.lat) / (C.lon - B.lon);
            T a_10 = B.lat - a_11 * B.lon;
            T b_1 = (C.lat - A.lat) / (C.lon - A.lon);
            T b_0 = A.lat - b_1 * A.lon;
            
            // test if lines are constructed correctly
            //assertlat(fequals(A.lon * a_01 + a_00, A.lat));
            //assertlat(fequals(B.lon * a_01 + a_00, B.lat));
            
            assertlat(fequals(A.lon * b_1 + b_0, A.lat));
            assertlat(fequals(C.lon * b_1 + b_0, C.lat));
            
            assertlat(fequals(B.lon * a_11 + a_10, B.lat));
            assertlat(fequals(C.lon * a_11 + a_10, C.lat));
            
            if(a_11 != 0)intA = 1.0 / a_11 * (cos(a_11 * B.lon + a_10) - cos(a_11 * C.lon + a_10));
            else intA = sin(a_10) * (C.lon - B.lon);
            
            if(b_1 != 0)intB = 1.0 / b_1 * (cos(b_1 * B.lon + b_0) - cos(b_1 * C.lon + b_0));
            else intB = sin(b_0) * (C.lon - B.lon);
            
            return abs(intA + intB) * R * R;
            
        }
        //                                               /|
        // last two lon equal, meaning triangle of form /_| e.g
        else if(B.lon == C.lon) {
            // depending if a_01 = 0 or b_1 = 0, form integrals
            T intA = 0;
            T intB = 0;
            
            T a_01 = (B.lat - A.lat) / (B.lon - A.lon);
            T a_00 = A.lat - a_01 * A.lon;
            T b_1 = (C.lat - A.lat) / (C.lon - A.lon);
            T b_0 = A.lat - b_1 * A.lon;
            
            // test if lines are constructed correctly
            assertlat(fequals(A.lon * a_01 + a_00, A.lat));
            assertlat(fequals(B.lon * a_01 + a_00, B.lat));
            
            assertlat(fequals(A.lon * b_1 + b_0, A.lat));
            assertlat(fequals(C.lon * b_1 + b_0, C.lat));
            
            // assertlat(fequals(B.lon * a_11 + a_10, B.lat));
            // assertlat(fequals(C.lon * a_11 + a_10, C.lat));
            
            if(a_01 != 0)intA = 1.0 / a_01 * (cos(a_01 * A.lon + a_00) - cos(a_01 * B.lon + a_00));
            else intA = sin(a_00) * (B.lon - A.lon);
            
            if(b_1 != 0)intB = 1.0 / b_1 * (cos(b_1 * A.lon + b_0) - cos(b_1 * C.lon + b_0));
            else intB = sin(b_0) * (C.lon - A.lon);
            
            return abs(intA + intB) * R * R;
            
        }
        // regular case, triangle is not aligned at any side to axis!
        else {
            
            //#error HERE CALCULATIONS ARE WRONG!
            
            
            // depending if a_01=0, a_11 = 0 or b_1 = 0, form integrals
            T intA1 = 0;
            T intA2 = 0;
            T intB = 0;
            
            T a_01 = (B.lat - A.lat) / (B.lon - A.lon);
            T a_00 = A.lat - a_01 * A.lon;
            T a_11 = (C.lat - B.lat) / (C.lon - B.lon);
            T a_10 = B.lat - a_11 * B.lon;
            T b_1 = (C.lat - A.lat) / (C.lon - A.lon);
            T b_0 = A.lat - b_1 * A.lon;
            
            // test if lines are constructed correctly
            assertlat(fequals(A.lon * a_01 + a_00, A.lat));
            assertlat(fequals(B.lon * a_01 + a_00, B.lat));
            
            assertlat(fequals(A.lon * b_1 + b_0, A.lat));
            assertlat(fequals(C.lon * b_1 + b_0, C.lat));
            
            assertlat(fequals(B.lon * a_11 + a_10, B.lat));
            assertlat(fequals(C.lon * a_11 + a_10, C.lat));
            
            if(a_11 != 0)intA1 = 1.0 / a_11 * (cos(a_11 * B.lon + a_10) - cos(a_11 * C.lon + a_10));
            else intA1 = sin(a_10) * (C.lon - B.lon);
            
            if(b_1 != 0)intB = 1.0 / b_1 * (cos(b_1 * C.lon + b_0) - cos(b_1 * A.lon + b_0));
            else intB = sin(b_0) * (C.lon - A.lon);
            
            if(a_01 != 0)intA2 = 1.0 / a_01 * (cos(a_01 * A.lon + a_00) - cos(a_01 * B.lon + a_00));
            else intA2 = sin(a_00) * (B.lon - A.lon);
            
            return abs(intA1 + intA2 + intB) * R * R;
            
        }
    }
    
    /// helper function which uses lltriangle to decompose a (convex!) polygon into several small triangles
    template<typename T> static T llPolygonArea(Vector3<T> *points, const unsigned int count, const double R = 1.0) {
        // make easy triangulation
        
        if(count == 0)return 0.0;
        
        // quicker round for 3 points
        if(count == 3)
        {
            Vector2<T> tA = points[0] / points[0].z;
            Vector2<T> tB = points[1] / points[1].z;
            Vector2<T> tC = points[2] / points[2].z;
            sortABClon(tA, tB, tC);
            return llTriangle(tA, tB, tC);
        }
        
        // first step is to sort points
        // therefore go through points and calc their angle in relation to y axes! (because it can be easily calculated with tan!)
        Vector3<T> center = Vector3<T>(0);
        for(int i = 0; i < count; i++) {
            if(points[i].z == 0)points[i].z = 1;
            points[i] = points[i] / points[i].z;
            center += points[i];
        }
        center = center / (double)count;
        center.z = 1;
        
        // dehomogenize
        //center = center / center.z;
        
        // calculate angles
        double *angles = new double[count];
        Vector3<T> m = cross(center, points[0]);
        Vector3<T> M = cross(Vector3<T>(0, 0, 1), m); // inf point!
        M = vector3Normalize(M);
        for(int i = 1; i < count; i++) {
            
            Vector3<T> l = cross(center, points[i]);
            
            // to retrieve angle, get inf points!
            Vector3<T> L = cross(Vector3<T>(0, 0, 1), l);
            
            // normalize both vectors
            L = vector3Normalize(L);
            
            // set angle, based on dot product of normals!
            angles[i] = (L * M);
            
            assert(-PI <= angles[i] && angles[i] <= PI);
        }
        
        // sort with insertion sort
        for(int i = 1; i < count; i++) {
            for(int k = i; k > 0 && angles[k] < angles[k - 1]; k--) {
                swap(points[k], points[k - 1]);
                swap(angles[k], angles[k - 1]);
            }
        }
        
        int pointcount = count;
        
        T area = 0;
        Vector3<T> A = points[0];
        Vector3<T> B;
        Vector3<T> C;
        
        Vector2<T> tA;
        Vector2<T> tB;
        Vector2<T> tC;
        // now calc area based on triangulation
        for(int i = 1; i < pointcount - 1; i++) {
            B = points[i];
            C = points[i + 1];
            
            tA = A;
            tB = B;
            tC = C;
            sortABClon(tA, tB, tC);
            
            area += Metrics::llTriangle(tA, tB, tC, R);
        }
        
        SafeDeleteA(angles);
        
        return area;
    }
    
    
    /// function to implement the volume
    template<typename T> static T linearFrustrum(const Vector2<T> _A, const Vector2<T> _B, const Vector2<T> _C, Vector2<T> *_points, T *depths, const unsigned int _w, const unsigned int _h, const double R = 1.0) {
        using namespace std;
        
        // (0) reorient if necessary, all data into NEMO Grid
        int w = _w;
        int h = _h;
        int orientation = Metrics::orientationOfCoords(_points, w, h);
        
        cout<<"before reorientation (NEMO Style):"<<endl;
        printArrayNEMOStyle(_points, w, h);
        cout<<"in mem:"<<endl;
        printArray(_points, w, h);
        Metrics::reorientArrayToNEMO(orientation, _points, w, h, true);
        
        cout<<"after reorientation:"<<endl;
        printArrayNEMOStyle(_points, w, h);
        cout<<"in mem:"<<endl;
        printArray(_points, w, h);
        w = _w;
        h = _h;
        Metrics::reorientArrayToNEMO(orientation, depths, w, h);
        assert(Metrics::orientationOfCoords(_points, w, h) == OM_NEMO);
        
        // (1) first step is to get bounding box triangle in lon/lat coords
        Vector2<T> minABCll = _A;
        Vector2<T> maxABCll = _B;
        minABCll = vector2Min(minABCll, _B);
        minABCll = vector2Min(minABCll, _C);
        maxABCll = vector2Max(maxABCll, _B);
        maxABCll = vector2Max(maxABCll, _C);
        
        // (2) find bounding box indices!
        int si=-1,ei=-1,sj=-1,ej=-1;
        Metrics::findBoundingBoxIndices(_points, w, h, minABCll, maxABCll, si, ei, sj, ej);
        
        // (3) now construct the points boundingBox
        //
        Vector2<T> pointsMin = _points[si + sj * w];
        Vector2<T> pointsMax = _points[ei + ej * w];
        
        // (3.1) if special case is given, that triangle encapsulates point region(which does not make sense!), alter pointgrid bounding cell and message out!
        if(pointsMin.lon > minABCll.lon || pointsMax.lon < maxABCll.lon || pointsMin.lat > minABCll.lat || pointsMax.lat < maxABCll.lat) {
            pointsMin = vector2Min(pointsMin, minABCll);
            pointsMax = vector2Max(pointsMax, maxABCll);
            cout<<">> info: all data points lie in triangle, are you sure?"<<endl;
        }
        
        // assert this is a real bounding Box!
        assert(pointsMin.lon <= pointsMax.lon);
        assert(pointsMin.lat <= pointsMax.lat);
        
        // (4) now convert all lon/lat coordinates into (homogenized) lon/lat coords
        // first triangle coords
        Vector3<T> A = Vector3<T>(_A.x, _A.y, 1);
        Vector3<T> B = Vector3<T>(_B.x, _B.y, 1);
        Vector3<T> C = Vector3<T>(_C.x, _C.y, 1);
        
        // first step, convert all lon/lat in xyz coordinates!
        Vector3<T> *points = new Vector3<T>[w * h];
        assert(points);
        for(int i = 0; i < w * h; i++) {
            points[i] = Vector3<T>(_points[i].x, _points[i].y);
        }
        
        // the volume...
        T volume = 0;
        
        // dehomogenized points in XYZ
        Vector2<T> trianglePoints[3];
        trianglePoints[0] = A / A.z;
        trianglePoints[1] = B / B.z;
        trianglePoints[2] = C / C.z;
        
        Vector2<T> p[4]; // edge points in XYZ
        
        Vector3<T> *areaPoints = new Vector3<T>[BUFFER_SIZE];
        int areaPointsCount = 0;
        
        // go through all cells (start with indices obtained in previous section)
        // as we access cells by i, i-1, j, j-1, shift indices by 1!
        si += 1; sj += 1; ei += 1; ej += 1;
        // loop through all points in BB
        for(int j=sj; j < ej; j++)
            for(int i=si; i < ei; i++)
            {
                // bottom layer can't be integrated
                if(i == 0 || j == 0)continue;
                
                // get Edge points
                p[0] = points[(i - 1)  + j * w]; // top left
                p[1] = points[i        + j * w]; // top right
                p[2] = points[i        + (j - 1) * w]; // bottom right
                p[3] = points[(i - 1)  + (j - 1) * w]; // bottom left
                
                // set back to zero
                areaPointsCount = 0;
                
                // now find out which points of A, B, C lie in Cell and vice versa
                for(int k = 0; k < 3; k++) {
                    if(pointInCell(p, trianglePoints[k]))areaPoints[areaPointsCount++] = Vector3<T>(trianglePoints[k].x, trianglePoints[k].y, 1);
                }
                for(int k = 0; k < 4; k++) {
#warning not yet implemented!!!
                    //if(pointInTriangle(trianglePoints, 3, p[k]))areaPoints[areaPointsCount++] = Vector3<T>(p[k].x, p[k].y, 1);
                }
                
                // now the most interesting part! determine intersection points and add to array
                quadrilateralTriangleIntersection(trianglePoints, p, areaPoints + areaPointsCount, &areaPointsCount);
                
                // print out points
                for(int u = 0; u < areaPointsCount; u++)cout<<"p"<<u<<": "<<areaPoints[u].toString()<<endl;
                
                // if less than 3 points skip! (would mean points are on one line!)
                if(areaPointsCount < 3)continue;
                
                // get area of spherical polygon, defined by points in vector areaPoints!
                double area = llPolygonArea(areaPoints, areaPointsCount, R);
                volume += area * depths[i + j * w] * 0.001;
            }
        
        SafeDeleteA(areaPoints);
        SafeDeleteA(points);
        return volume;
    }
    
    /// function to implement the area metric
    template<typename T> static T linearArea(const Vector2<T> _A, const Vector2<T> _B, Vector2<T> *_points, T *depths, const unsigned int _w, const unsigned int _h, const double R = 1.0) {
        using namespace std;
        
        // (0) reorient if necessary, all data into NEMO Grid
        int w = _w;
        int h = _h;
        int orientation = Metrics::orientationOfCoords(_points, w, h);
        Metrics::reorientArrayToNEMO(orientation, _points, w, h);
        w = _w;
        h = _h;
        Metrics::reorientArrayToNEMO(orientation, depths, w, h);
        assert(Metrics::orientationOfCoords(_points, w, h) == OM_NEMO);
        
        // (1) first step is to get bounding box triangle in lon/lat coords
        Vector2<T> minABll = _A;
        Vector2<T> maxABll = _B;
        minABll = vector2Min(minABll, _B);
        maxABll = vector2Max(maxABll, _B);
        
        // (2) find bounding box indices!
        int si=-1,ei=-1,sj=-1,ej=-1;
        Metrics::findBoundingBoxIndices(_points, w, h, minABll, maxABll, si, ei, sj, ej);
        
        // (3) now construct the points boundingBox
        //
        Vector2<T> pointsMin = _points[si + sj * w];
        Vector2<T> pointsMax = _points[ei + ej * w];
        
        // assert this is a real bounding Box!
        assert(pointsMin.lon <= pointsMax.lon);
        assert(pointsMin.lat <= pointsMax.lat);
        
        // (4) now convert all lon/lat coordinates into (homogenized) lon/lat coords
        // first triangle coords
        Vector3<T> A = Vector3<T>(_A.x, _A.y, 1);
        Vector3<T> B = Vector3<T>(_B.x, _B.y, 1);
        
        // calc line between A & B
        Vector3<T> line = cross(A, B);
        
        // first step, convert all lon/lat in xyz coordinates!
        Vector3<T> *points = new Vector3<T>[w * h];
        assert(points);
        for(int i = 0; i < w * h; i++) {
            points[i] = Vector3<T>(_points[i].x, _points[i].y);
        }
        
        // the area...
        T area = 0;
        
        // dehomogenized points in XYZ
        Vector2<T> lPoints[2];
        lPoints[0] = A / A.z;
        lPoints[1] = B / B.z;
        
        Vector2<T> p[4]; // edge points in XYZ
        
        Vector3<T> *linePoints = new Vector3<T>[BUFFER_SIZE];
        int linePointsCount = 0;
        
        // go through all cells (start with indices obtained in previous section)
        // as we access cells by i, i-1, j, j-1, shift indices by 1!
        si += 1; sj += 1; ei += 1; ej += 1;
        // loop through all points in BB
        for(int j=sj; j < ej; j++)
            for(int i=si; i < ei; i++)
            {
                // bottom layer can't be integrated
                if(i == 0 || j == 0)continue;
                
                // get Edge points
                p[0] = points[(i - 1)  + j * w]; // top left
                p[1] = points[i        + j * w]; // top right
                p[2] = points[i        + (j - 1) * w]; // bottom right
                p[3] = points[(i - 1)  + (j - 1) * w]; // bottom left
                
                // set back to zero
                linePointsCount = 0;
                
                // now find out which points of A, B lie in Cell
                for(int k = 0; k < 2; k++) {
                    if(pointInCell(p, lPoints[k]))linePoints[linePointsCount++] = Vector3<T>(lPoints[k].x, lPoints[k].y, 1);
                }
                
                // now the most interesting part! determine intersection points and add to array
                quadrilateralLineIntersection(A, B, p, linePoints + linePointsCount, &linePointsCount);
                
                
                // have two intersection points been found?
                assert(linePointsCount <= 2);
                if(linePointsCount == 2) {
                    // note that depth is in m^2! distance should return something in km!
                    area += Metrics::llDistance(DHom(linePoints[0]), DHom(linePoints[1])) * depths[i + j * w] * 0.001;
                }
            }
        
        SafeDeleteA(linePoints);
        SafeDeleteA(points);
        return area;
    }
    
    /// helper to get derivative of spline
    /// destructive
    template<typename T> static void polyDerivative(T *poly, unsigned int n) {
        for(int i = 1; i < n; i++) {
            poly[i - 1] = poly[i] * i;
        }
        poly[n - 1] = 0;
    }
    
    
    /// helper which gives intersection point with a spline and a line
    /// splineA, splineB define the interval for which the spline is defined
    /// A, B define a line
    template<typename T> static bool    splineLineIntersectionLon(const Vector2<T> A, const Vector2<T> B, const T *spline, T splineALon, T splineBLon, Vector3<T> *out, unsigned int *count) {
        
        bool intersection = false;
        
        // order splineAlon, splineBlon
        if(splineALon > splineBLon)swap(splineBLon, splineALon);
        assert(splineALon <= splineBLon);
        
        // first of all decide two cases!
        // y = ax +b or x = d
        if(fequals(A.lon, B.lon)) {
            // easy case, x = d -- where d == A.lon!
            
            // if A.lon is in interval splineA/splineB intersection points can be directly calculated by evaluating the spline!
            if(A.lon >= splineALon && A.lon <= splineBLon) {
                // lies in interval! --> valid intersection point
                Vector3<T> IP = Vector3<T>(A.lon, evalPolynomial(A.lon, spline, 4), 1);
                // check now if also on valid segment of line A <-> B
                if((IP.lat >= A.lat && IP.lat <= B.lat) || (IP.lat <= A.lat && IP.lat >= B.lat)) { // this check is doubled, cause ordering of A,B should not matter!
                    *(out) = IP;
                    (*count)++;
                    intersection = true;
                }
            }
            
        }
        else {
            // g(x) = ax + b
            // construct a & b
            T a = (B.lat - A.lat) / (B.lon - A.lon);
            T b = A.lat - a * A.lon;
            
            // assert construction was ok
            assert(fequals(a * A.lon + b, A.lat));
            assert(fequals(a * B.lon + b, B.lat));
            
            // now find cubic roots of polynomial s(x)-g(x) (where s(x) is the spline's cubic polynomial!)
            T s[4] = {spline[0], spline[1], spline[2], spline[3]};
            s[0] -= b;
            s[1] -= a;
            
            // to prevent numerical errors, set values near to zero, to zero
            for(int i = 0; i < 4; i++)if(fequals(s[i], 0))s[i] = 0;
            
            T roots[3]; // there can be 3 roots!
            int numRoots = 0;
            cubicRealRoots(s[3], s[2], s[1], s[0], roots, &numRoots);
            
            int pos = 0;
            // test if root lies in interval!
            for(int i = 0; i < numRoots; i++) {
                // use epsilon for special case allowance!
                if(roots[i] >= splineALon - EPSILON && roots[i] <= splineBLon + EPSILON) {
                    // lies in interval! --> valid intersection point
                    
                    // test now, if points lies also in interval defined by points A & B!
                    if((roots[i] >= A.lon - EPSILON && roots[i] <= B.lon + EPSILON) || (roots[i] <= A.lon + EPSILON && roots[i] >= B.lon - EPSILON)) {
                        Vector3<T> IP = Vector3<T>(roots[i], a * roots[i] + b, 1);
                        *(out + pos) = IP;
                        pos++;
                        (*count)++;
                        intersection = true;
                    }
                }
            }
            
        }
        return intersection;
    }
    
    /// the function is the same for lat, but we flip it right!
    template<typename T> static bool    splineLineIntersectionLat(const Vector2<T> A, const Vector2<T> B, const T *spline, const T splineALat, const T splineBLat, Vector3<T> *out, unsigned int *count) {
        
        // save first, because lon/lats have to flipped back after computation
        Vector3<T>      points[3]; // no more than 3 intersection points!
        unsigned int    pointscount = 0;
        bool res = splineLineIntersectionLon(Vector2<T>(A.y, A.x), Vector2<T>(B.y, B.x), spline, splineALat, splineBLat, points, &pointscount);
        
        // copy & flip coords, if no error occured
        if(res) {
            for(int i = 0; i < pointscount; i++) {
                // flip coordinates
                *(out + i) = Vector3<T>(points[i].y, points[i].x, 1);
                (*count)++;
            }
        }
        
        return res;
    }
    
    /// function to find bounding box indices
    template<typename T> static void findBoundingBoxIndices(Vector2<T> * _points, const int w, const int h, const Vector2<T> min, const Vector2<T> max, int& si, int & ei, int& sj, int & ej) {
        // (2) second step is to find indices of bounding box of points
        si = w - 1;
        ei = 0;
        sj = h - 1;
        ej = 0;
        
        // (2.1) get start / end indices in i direction
        // general idea is here to get the first cell from front/back,
        // where the nearer edge is inside the line's BoundingBox
        // the +1/-1 test the nearer edge
        // after that an assert loop is run, because this has to be secured later
        
        // get start indices for lon
        for(int j = 0; j < h; j++)
        {
            int cursi = 0;
            int curei = w - 1;
            
            // find current si
            while(_points[(cursi + 1) + j * w].lon < min.lon)cursi++;
            // find current ei
            while(_points[(curei - 1) + j * w].lon > max.lon)curei--;
            
            // global min/max along all columns!
            si = ::min(si, cursi);
            ei = ::max(ei, curei);
        }
#warning write better in area functions a checker if line/triangle or whatever is outside of datapoints!
#warning assert makes only sense if si != 0, ...
        // assert everything is ok!
        for(int j = 0; j < h; j++) {
            //  assert(_points[si + j * w].lon <= min.lon);
            //assert(_points[ei + j * w].lon >= max.lon);
        }
        
        // now for lat
        for(int i = 0; i < w; i++)
        {
            int cursj = 0;
            int curej = h - 1;
            
            // find current si
            while(_points[i + (cursj + 1) * w].lat < min.lat)cursj++;
            // find current ei
            while(_points[i + (curej - 1) * w].lat > max.lat)curej--;
            
            // global min/max along all columns!
            sj = ::min(sj, cursj);
            ej = ::max(ej, curej);
        }
        
#warning write better in area functions a checker if line/triangle or whatever is outside of datapoints!
#warning assert makes only sense if sj != 0, ...
        // assert everything is ok!
        for(int i = 0; i < w; i++) {
            //assert(_points[i + sj * w].lat <= min.lat);
            //assert(_points[i + ej * w].lat >= max.lat);
        }
    }
    
    /// helper function to solve the special case when there are two/three intersections per side
    template<typename T> static void removePointsToPreventDoubleCounting(const Vector2<T> *cellpoints, Vector3<T> *points, unsigned int start, unsigned int& count) {
        
        // points get removed if inside polygon!
        int pos = start;
        while(pos < count) {
            // point inside polygon?
            if(pointInPoly(cellpoints, 4, points[pos])) {
                // delete point from list
                std::swap(points[pos], points[--count]);
            }
            else {
                // inc pointer
                pos++;
            }
        }
    }
    
    /// helper function to solve integral between cubic polynomial of spline and line constructed from 2 points
    /// a, b are the limits
    template<typename T> static T splineLineIntegral(const T a, const T b, const T *spline, const Vector2<T> A, const Vector2<T> B) {
        
        // not a necessary assertion...
        //assert(a < b);
        
        T aa = a * a;
        T bb = b * b;
        T aaa = a * aa;
        T bbb = b * bb;
        T aaaa = a * aaa;
        T bbbb = b * bbb;
        
        assert(!fequals(A.lon, B.lon));
        
        T m = (B.lat - A.lat) / (B.lon - A.lon);
        T c = A.lat - m * A.lon;
        
        // assert construction was ok
        assert(fequals(m * A.lon + c, A.lat));
        assert(fequals(m * B.lon + c, B.lat));
        
        // now find cubic roots of polynomial s(x)-g(x) (where s(x) is the spline's cubic polynomial!)
        T s[4] = {spline[0], spline[1], spline[2], spline[3]};
        s[0] -= c;
        s[1] -= m;
        
        return (0.25 * s[3] * bbbb + 1/3.0 * s[2] * bbb + 0.5 * s[1] * bb + s[0] * b) -
        (0.25 * s[3] * aaaa + 1/3.0 * s[2] * aaa + 0.5 * s[1] * aa + s[0] * a);
    }
    
    /// function to return all line segments, with value for splineArea intersection
    template<typename T> static std::vector<SLineSegment<T> > splineSegmentsInsideCell(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& topleft, const Vector2<T>& topright, const Vector2<T>& bottomright, const Vector2<T>& bottomleft, T *xTopSpline, T *xBottomSpline, T *yLeftSpline, T *yRightSpline) {
        using namespace std;
        vector<SLineSegment<T> > res;
        
        // some asserts
        assert(xTopSpline);
        assert(xBottomSpline);
        assert(yLeftSpline);
        assert(yRightSpline);
        
        // intersect with splines
        // all in all,
        Vector3<T> intersectionPoints[4 * 3]; // 0-2 top, 3-5 bottom, 6-8 left, 9-11 right
        unsigned int        intersectionPointsCount[4] = {0, 0, 0, 0}; // 0 top, 1 bottom, 2 left, 3 right
        
        Metrics::splineLineIntersectionLon(A, B, xTopSpline, topleft.x, topright.x,
                                           intersectionPoints + SIDE_TOP * 3, &intersectionPointsCount[SIDE_TOP]);
        Metrics::splineLineIntersectionLon(A, B, xBottomSpline, bottomleft.x, bottomright.x,
                                           intersectionPoints + SIDE_BOTTOM * 3, &intersectionPointsCount[SIDE_BOTTOM]);
        Metrics::splineLineIntersectionLat(A, B, yLeftSpline, topleft.y, bottomleft.y,
                                           intersectionPoints + SIDE_LEFT * 3, &intersectionPointsCount[SIDE_LEFT]);
        Metrics::splineLineIntersectionLat(A, B, yRightSpline, topright.y, bottomright.y,
                                           intersectionPoints + SIDE_RIGHT * 3, &intersectionPointsCount[SIDE_RIGHT]);
        
        // now add all those points to array with sides!
        int intersectionPointsWithSidesCount = 0;
        SIntersectionPoint<T> intersectionPointsWithSides[14];
        for(int k = 0; k < 4; k++) {
            for(int l = 0; l < intersectionPointsCount[k]; l++) {
                intersectionPointsWithSides[intersectionPointsWithSidesCount].point = DHom(intersectionPoints[l + k * 3]);
                intersectionPointsWithSides[intersectionPointsWithSidesCount].side = (CellSide)k;
                intersectionPointsWithSidesCount++;
            }
        }
        
        // now find out which points of A, B lie in Cell
            if(Metrics::llPointInSplineCell(A, yLeftSpline, yRightSpline, xTopSpline, xBottomSpline)) {
                intersectionPointsWithSides[intersectionPointsWithSidesCount].point = A;
                intersectionPointsWithSides[intersectionPointsWithSidesCount].side = SIDE_UNKNOWN;
                intersectionPointsWithSidesCount++;
            }
        if(Metrics::llPointInSplineCell(B, yLeftSpline, yRightSpline, xTopSpline, xBottomSpline)) {
            intersectionPointsWithSides[intersectionPointsWithSidesCount].point = B;
            intersectionPointsWithSides[intersectionPointsWithSidesCount].side = SIDE_UNKNOWN;
            intersectionPointsWithSidesCount++;
        }
        //
        // sort intersection points with sides array after lon(if lons are equal sort after lat!)
        //qsort((void*)intersectionPointsWithSides, intersectionPointsWithSidesCount, sizeof(SIntersectionPoint<T>), compIntersectionPoints<T>);
        
        
        // calculate distance in cell!
        
        T distance = 0;
        SLineSegment<T> ls;
        // go through all pairs, and decide if their distance counts!
        for(int k = 0; k < intersectionPointsWithSidesCount - 1; k++) {
            // pairs
            // if the pair lies on different sides, it counts for sure!
            // (follows from Jordan's curve theorem!)
            if(intersectionPointsWithSides[k].side != intersectionPointsWithSides[k + 1].side) {
                //distance += Metrics::llDistance(intersectionPointsWithSides[k].point, intersectionPointsWithSides[k + 1].point);
                ls.X = intersectionPointsWithSides[k].point; ls.Y = intersectionPointsWithSides[k + 1].point;
                res.push_back(ls);
            }
            else {
                // more difficult case!, both points are on same side --> decide with help of integral if connecting line is outside or inside cell!
                // two points, decide now with help of integral if distance should be added or not!
                Vector2<T> pointA = intersectionPointsWithSides[k].point;
                Vector2<T> pointB = intersectionPointsWithSides[k + 1].point;
                CellSide side = intersectionPointsWithSides[k].side;
                
                // now decide depending on sign of integral, if line between points lie in Cell or not!
                // x-based sides
                if(side == SIDE_TOP) {
                    // top spline
                    T intA = pointA.lon < pointB.lon ? pointA.lon : pointB.lon;
                    T intB = pointA.lon < pointB.lon ? pointB.lon : pointA.lon;
                    T integral = Metrics::splineLineIntegral(intA, intB, xTopSpline, pointA, pointB);
                    
                    // if integral > 0 ==> add to distance!
                    if(integral > 0) { //distance += Metrics::llDistance(pointA, pointB);
                        ls.X = pointA; ls.Y = pointB;
                        res.push_back(ls);
                    }
                    
                }
                if(side == SIDE_BOTTOM) {
                    // bottom spline
                    T intA = pointA.lon < pointB.lon ? pointA.lon : pointB.lon;
                    T intB = pointA.lon < pointB.lon ? pointB.lon : pointA.lon;
                    T integral = Metrics::splineLineIntegral(intA, intB, xBottomSpline, pointA, pointB);
                    
                    // if integral < 0 ==> add to distance!
                    if(integral < 0) { //distance += Metrics::llDistance(pointA, pointB);
                        ls.X = pointA; ls.Y = pointB;
                        res.push_back(ls);
                    }
                }
                // y based sides
                if(side == SIDE_LEFT) {
                    // left spline
                    T intA = pointA.lat < pointB.lat ? pointA.lat : pointB.lat;
                    T intB = pointA.lat < pointB.lat ? pointB.lat : pointA.lat;
                    T integral = Metrics::splineLineIntegral(intA, intB, yLeftSpline, swapXY(pointA), swapXY(pointB));
                    
                    // if integral > 0 ==> add to distance!
                    if(integral > 0) { //distance += Metrics::llDistance(pointA, pointB);
                        ls.X = pointA; ls.Y = pointB;
                        res.push_back(ls);
                    }
                    
                }
                if(side == SIDE_RIGHT) {
                    // right spline
                    T intA = pointA.lat < pointB.lat ? pointA.lat : pointB.lat;
                    T intB = pointA.lat < pointB.lat ? pointB.lat : pointA.lat;
                    T integral = Metrics::splineLineIntegral(intA, intB, yRightSpline, swapXY(pointA), swapXY(pointB));
                    
                    // if integral < 0 ==> add to distance!
                    if(integral < 0) { //distance += Metrics::llDistance(pointA, pointB);
                        ls.X = pointA; ls.Y = pointB;
                        res.push_back(ls);
                    }
                }
                
                // A & B directly in one cell
                if(side == SIDE_UNKNOWN) {
                    //distance += Metrics::llDistance(pointA, pointB);
                    ls.X = pointA; ls.Y = pointB;
                    res.push_back(ls);
                }
                
            }
        }
        
        return res;
    }
    
    /// function to implement the area metric using cubic splines as approximation for cells.
    template<typename T> static T splineArea(const Vector2<T> _A, const Vector2<T> _B, Vector2<T> *_points, T *depths, const unsigned int _w, const unsigned int _h, const double R = 1.0) {
        using namespace std;
        
        // output filestream
        ofstream file;
        file.open("/Users/leonhardspiegelberg/Documents/MATLAB/testgrid.m");
        file<<"% c++ generated matlab file, to plot calculations\n";
        file<<"fig = figure();\n";
        
        // output filestream for tikzimage
        ofstream tfile;
        tfile.open("/Users/leonhardspiegelberg/Documents/UniversityofAlberta/Bachelorarbeit/thesis/img/gridtex.tex");
        tfile<<"%% tex file for a tikz image"<<endl;
        tfile<<"\\begin{tikzpicture}[scale=0.1]"<<endl;
        tfile.precision(5);
        tfile.setf(ios::fixed);
        
        // output file to write out all line segments, with lengths
        ofstream cfile;
        cfile.open("/Users/leonhardspiegelberg/Documents/MATLAB/segments.txt");
        // print format info, first
        cfile<<"# file storing all line segments in cells"<<endl;
        cfile<<"# only cells containing at least one line segment are printed out with format:"<<endl;
        cfile<<"# begin cell(i, j):"<<endl;
        cfile<<"# segmentlengthsum"<<endl;
        cfile<<"# firstpointofsegment secondpointofsegment segmentlength"<<endl;
        cfile<<"# ..."<<endl;
        cfile<<"# end cell(i, h):"<<endl;
        
        // red color table
        double3 ctable[] = {double3(255, 48, 48),
                            double3(238, 44, 44),
                            double3(205, 38, 38),
                            double3(139, 26, 26),
                            double3(255, 64, 64),
                            double3(238, 59, 59),
                            double3(205, 51, 51),
                            double3(139, 35, 35),
                            double3(255, 140, 105),
                            double3(139, 76, 57),
                            double3(205, 0, 0),
                            double3(139, 0, 0)};
        int ctable_count = 12;
        for(int i = 0; i < ctable_count; i++)ctable[i] *= 1/255.0;
        
        // (0) reorient if necessary, all data into NEMO Grid
        int w = _w;
        int h = _h;
        int orientation = Metrics::orientationOfCoords(_points, w, h);
        Metrics::reorientArrayToNEMO(orientation, _points, w, h);
        w = _w;
        h = _h;
        Metrics::reorientArrayToNEMO(orientation, depths, w, h);
        assert(Metrics::orientationOfCoords(_points, w, h) == OM_NEMO);
        
        // (1) first step is to get bounding box triangle in lon/lat coords
        Vector2<T> minABll = _A;
        Vector2<T> maxABll = _B;
        minABll = vector2Min(minABll, _B);
        maxABll = vector2Max(maxABll, _B);
        
        // (2) find bounding box indices!
        int si=-1,ei=-1,sj=-1,ej=-1;
        Metrics::findBoundingBoxIndices(_points, w, h, minABll, maxABll, si, ei, sj, ej);
        
        // (3) now construct the points boundingBox
        //
        Vector2<T> pointsMin = _points[si + sj * w];
        Vector2<T> pointsMax = _points[ei + ej * w];
        
        // assert this is a real bounding Box!
        assert(pointsMin.lon <= pointsMax.lon);
        assert(pointsMin.lat <= pointsMax.lat);
        
        // (4) now convert all lon/lat coordinates into (homogenized) lon/lat coords
        // first triangle coords
        Vector3<T> A = Vector3<T>(_A.x, _A.y, 1);
        Vector3<T> B = Vector3<T>(_B.x, _B.y, 1);
        
        // calc line between A & B
        Vector3<T> line = cross(A, B);
        
        // first step, convert all lon/lat in xyz coordinates!
        //Vector3<T> *points = new Vector3<T>[w * h];
        // NEW! use Grid!
        TArray<Vector3<T> > points(w, h);
        //  points.setMode(GM_MIRRORX | GM_CLAMPY);
        points.setMode(GM_CLAMP);
        for(int i = 0; i < w * h; i++) {
            points.set(i,Vector3<T>(_points[i].x, _points[i].y));
        }
        
        // the area...
        T area = 0;
        
        // dehomogenized points in XYZ
        Vector2<T> lPoints[2];
        lPoints[0] = A / A.z;
        lPoints[1] = B / B.z;
        
        Vector2<T> p[4]; // edge points in XYZ
        
        Vector3<T> *linePoints = new Vector3<T>[BUFFER_SIZE];
        unsigned int linePointsCount = 0;
        
        // go through all cells (start with indices obtained in previous section)
        // as we access cells by i, i-1, j, j-1, shift indices by 1!
        si += 1; sj += 1; ei += 1; ej += 1;
        
        // ensure that there are enough points!
        assert(si < ei);
        assert(sj < ej);
        
        // line
        file<<"plot(["<<A.x<<" "<<B.x<<"], ["<<A.y<<","<<B.y<<"]);\nhold on;\n%ending c++ generated code\n"<<flush;
        
#warning add here a better check for boundaries!
        
        // buffer for all intersection points
        SIntersectionPoint<T> intersectionPointsWithSides[14]; // maximum 14 needed!
        int intersectionPointsWithSidesCount = 0;
        
        // loop through all points in BB
        for(int j=sj; j < ej; j++)
            for(int i=si; i < ei; i++)
            {
                // bottom layer can't be integrated
                if(i == 0 || j == 0)continue;
                
                // get Edge points
                p[0] = points.get(i - 1, j    ); // top left
                p[1] = points.get(i    , j    ); // top right
                p[2] = points.get(i    , j - 1); // bottom right
                p[3] = points.get(i - 1, j - 1); // bottom left
                
                
                // construct the four splines describing the cell
                T yLeftSpline[4]; // function lat -> lon
                T yRightSpline[4]; // function lat -> lon
                T xTopSpline[4]; // function lon->lat
                T xBottomSpline[4]; // function lon->lat
                Vector2<T> s0, s1, s2, s3;
                
                // !!! important, new selector for grid points should be implemented!
                
                
                /// normal orientation(attention! order is important!!!)
                // left side
                s0 = points.get(i - 1, j - 2);
                s1 = points.get(i - 1, j - 1);
                s2 = points.get(i - 1, j + 0);
                s3 = points.get(i - 1, j + 1);
                cubicSplineCoefficients(s0.y, s1.y, s2.y, s3.y, s0.x, s1.x, s2.x, s3.x, yLeftSpline);
                
                // right side
                s0 = points.get(i + 0, j - 2);
                s1 = points.get(i + 0, j - 1);
                s2 = points.get(i + 0, j + 0);
                s3 = points.get(i + 0, j + 1);
                cubicSplineCoefficients(s0.y, s1.y, s2.y, s3.y, s0.x, s1.x, s2.x, s3.x, yRightSpline);
                
                // top side
                s0 = points.get(i - 2, j);
                s1 = points.get(i - 1, j);
                s2 = points.get(i - 0, j);
                s3 = points.get(i + 1, j);
                cubicSplineCoefficients(s0.x, s1.x, s2.x, s3.x, s0.y, s1.y, s2.y, s3.y, xTopSpline);
                
                // bottom side
                s0 = points.get(i - 2, j - 1);
                s1 = points.get(i - 1, j - 1);
                s2 = points.get(i - 0, j - 1);
                s3 = points.get(i + 1, j - 1);
                cubicSplineCoefficients(s0.x, s1.x, s2.x, s3.x, s0.y, s1.y, s2.y, s3.y, xBottomSpline);
                
                
                // determine all line segments which intersect with cell sides
                vector<SLineSegment<T> > ls = Metrics::splineSegmentsInsideCell(_A, _B, p[0], p[1], p[2], p[3],
                                                                      xTopSpline, xBottomSpline, yLeftSpline, yRightSpline);
                
                
                // calc distance
                T distance = 0;
                if(!ls.empty())
                    for(typename vector<SLineSegment<T> >::const_iterator it = ls.begin(); it != ls.end(); ++it) {
                        distance += Metrics::llDistance(it->X, it->Y);
                        
                        // print MATLAB code
                        double3 col = ctable[rand() % ctable_count];
                        double2 pointA = it->X;
                        double2 pointB = it->Y;
                        file<<"scatter("<<pointA.x<<","<<pointA.y<<", 15, 'o', 'fill', 'MarkerFaceColor',["<<col.x<<" "<<col.y<<" "<<col.z<<"], 'MarkerEdgeColor', 'none'); hold on;\n";
                        file<<"scatter("<<pointB.x<<","<<pointB.y<<", 15, 'o', 'fill', 'MarkerFaceColor',["<<col.x<<" "<<col.y<<" "<<col.z<<"], 'MarkerEdgeColor', 'none'); hold on;\n";
                        //file<<"plot(["<<pointA.x<<" "<<pointB.x<<"],["<<pointA.y<<" "<<pointB.y<<"], 'Color',["<<col.x<<" "<<col.y<<" "<<col.z<<"], 'LineWidth', 1); hold on; \n";
                    }
                
                // print out all data of segments
                if(!ls.empty()) {
                    // print out cell with index i,j)
#warning print out here original indices!
                    cfile<<"begin cell("<<i<<","<<j<<"):\n";
                    cfile<<distance<<endl;
                    for(typename vector<SLineSegment<T> >::const_iterator it = ls.begin(); it != ls.end(); ++it) {
                        cfile<<"("<<it->X.x<<","<<it->X.y<<") \t("<<it->Y.x<<","<<it->Y.y<<") \t"<<Metrics::llDistance(it->X, it->Y)<<endl;
                    }
                    cfile<<"end cell("<<i<<","<<j<<")\n";
                }
                    
                
                // print spline codes to matlab file
                file<<printOutSplineCode(xTopSpline, p[0].x, p[1].x, 'x');
                file<<printOutSplineCode(xBottomSpline, p[3].x, p[2].x, 'x');
                file<<printOutSplineCode(yLeftSpline, p[0].y, p[3].y, 'y');
                file<<printOutSplineCode(yRightSpline, p[1].y, p[2].y, 'y');
                
                // print splines into tikz file
                double globalscale = 0.1;
                tfile<<printOutTikzSplineCode(xTopSpline, p[0].x, p[1].x, 'x', globalscale);
                tfile<<printOutTikzSplineCode(xBottomSpline, p[3].x, p[2].x, 'x', globalscale);
                tfile<<printOutTikzSplineCode(yLeftSpline, p[0].y, p[3].y, 'y', globalscale);
                tfile<<printOutTikzSplineCode(yRightSpline, p[1].y, p[2].y, 'y', globalscale);
                
                // note that depth is in m^2! distance should return something in km!
                area += distance * depths[i + j * w] * 0.001;
                
                // flush stream
                file<<flush;
            }
        
        
        //tfile<<"\draw[scale=0.5,domain=-3.141:3.141,smooth,variable=\t] plot ({\t *\t},{\t});"<<endl;
        tfile<<"\\end{tikzpicture}"<<endl;
        tfile.close();
        
        file<<"set(gca, 'Visible', 'off')"<<endl;
        file<<"set(gcf, 'Color', 'white')"<<endl;
        //close file
        file.close();
        cfile.close();
        
        SafeDeleteA(linePoints);
        //SafeDeleteA(points);
        return area;
    }
    
    template<typename T> static std::string printOutTikzSplineCode(T *a, T start, T end, char mode, T globalscale = 1.0) {
         using namespace std;
        stringstream strstr;
        
        strstr.precision(5);
        strstr.setf(ios::fixed);
        
        // flip start/end if necessary
        if(start > end)::swap(start, end);
        start *= globalscale;
        end *= globalscale;
        
        strstr<<"\\draw [thin] plot [parametric,domain="<<start<<":"<<end;
        strstr<<",samples=50] function {";
        if(mode == 'x') {
            //t,t * t * t * -0.00001+ t * t * -0.00083+ t * -0.00066+9.99967
            strstr<<"t"<<",";
            strstr<<"t * t * t *"<<globalscale * a[3]<<" + ";
            strstr<<"t * t *"<<globalscale * a[2]<<" + ";
            strstr<<"t *"<<globalscale * a[1]<<" + ";
            strstr<<globalscale * a[0];
        }else if(mode == 'y') {
            strstr<<"t * t * t *"<<globalscale * a[3]<<" + ";
            strstr<<"t * t *"<<globalscale * a[2]<<" + ";
            strstr<<"t *"<<globalscale * a[1]<<" + ";
            strstr<<globalscale * a[0];
            strstr<<",t";
        }
        else {
            strstr<<"error!"<<endl;
        }
        strstr<<"};"<<endl;
        return strstr.str();
    }
    
    template<typename T> static std::string printOutSplineCode(T *a, T start, T end, char mode) {
        using namespace std;
        T h = 0.1;
        stringstream xcords;
        stringstream ycords;
        
        //flip if necessary
        if(start > end)std::swap(start, end);
        
        for(int i = 0; i < (end - start) / h; i++) {
            // eval spline
            T t = start + i * h;
            T val = evalPolynomial(t, a, 4);
            
            if(mode == 'x') {
                xcords<<t;
                ycords<<val;
                xcords<<" ";
                ycords<<" ";
            }
            else if(mode == 'y') {
                xcords<<val;
                ycords<<t;
                xcords<<" ";
                ycords<<" ";
            }
            else cout<<">> error: please specify 'x' or 'y' as mode"<<endl;
        }
        
        // last point
        // eval spline
        T t = end;
        T val = evalPolynomial(t, a, 4);
        
        if(mode == 'x') {
            xcords<<t;
            ycords<<val;
            xcords<<" ";
            ycords<<" ";
        }
        else if(mode == 'y') {
            xcords<<val;
            ycords<<t;
            xcords<<" ";
            ycords<<" ";
        }
        else cout<<">> error: please specify 'x' or 'y' as mode"<<endl;
        
        
        // print final code
        stringstream finalstr;
        finalstr<<"plot(["<<xcords.str()<<"], ["<<ycords.str()<<"], 'Color', [0.8 0.8 0.8]); hold on;\n";
        return finalstr.str();
    }
};
#endif /* defined(__nemoExtractor__Metrics__) */
