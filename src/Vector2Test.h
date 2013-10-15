//
//  Vector2Test.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 24.05.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef nemoExtractor_Vector2Test_h
#define nemoExtractor_Vector2Test_h

//CPP Unit
#include <cppunit/TestCase.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "MathFunctions.h"


// see http://cppunit.sourceforge.net/doc/lastest/cppunit_cookbook.html for details

/// Vector2 class tests
class Vector2Test : public CppUnit::TestFixture
{
private:
	//fast setup macros
	CPPUNIT_TEST_SUITE(Vector2Test);
	CPPUNIT_TEST(testOperators);
    CPPUNIT_TEST(testNearestNeighbour);
    CPPUNIT_TEST(testBilinearInterpolation);
    CPPUNIT_TEST(testPointInPoly);
    CPPUNIT_TEST(testLLDistance);
    CPPUNIT_TEST(testLLTriangleAndPoly);
    CPPUNIT_TEST(testPartRatio);
    CPPUNIT_TEST(testQuadIntersection);
    CPPUNIT_TEST(testSort3El);
    CPPUNIT_TEST(testIndexTransformations);
	CPPUNIT_TEST_SUITE_END();
    
    
public:
    
    /// test +, - , *, /
    void testOperators() {
        Vector2<double> v = Vector2<double>(2.0, 5.0);
        Vector2<double> v2 = Vector2<double>(2.0, 5.0);
        
        CPPUNIT_ASSERT(v == v2);
        CPPUNIT_ASSERT(v.equals(v2));
        v2.x = 2.0 + EPSILON / 2;
        
        CPPUNIT_ASSERT(v == v2);
        CPPUNIT_ASSERT(v.equals(v2));
        
        CPPUNIT_ASSERT(v + v2 == Vector2<double>(4.0, 10.0));
        CPPUNIT_ASSERT(v - v2 == Vector2<double>(0.0, 0.0));
        CPPUNIT_ASSERT(5.0 * v == Vector2<double>(10.0, 25.0));
        CPPUNIT_ASSERT(v * 5.0 == Vector2<double>(10.0, 25.0));
        CPPUNIT_ASSERT(v / 2.0 == Vector2<double>(1.0, 2.5));
        
        // length & dot product
        CPPUNIT_ASSERT(abs(v * v - 29.0) < EPSILON);
        v = Vector2<double>(4.0, 4.0);
        CPPUNIT_ASSERT(abs(v.length() - sqrt(4.0 * 4.0 + 4.0 * 4.0))<EPSILON);
    }
    
    /// test nearest neighbour Interpolation
    void testNearestNeighbour() {
        double a = 3.0;
        double b = 5.0;
        double c = 7.0;
        double d = 11.0;
        double x = 0.0;
        double y = 0.0;
        
        double2 A = double2(-2.0, 2.0);
        double2 B = double2(2.0, 2.0);
        double2 C = double2(-2.0, -2.0);
        double2 D = double2(2.0, -2.0);
        
        // perform four queries for each points
        x = -1.0;   y = 1.0;
        CPPUNIT_ASSERT(vector2NearestNeighbourEuclidean(A, B, C, D, a, b, c, d, x, y) == a);
        x = 1.0;   y = 1.0;
        CPPUNIT_ASSERT(vector2NearestNeighbourEuclidean(A, B, C, D, a, b, c, d, x, y) == b);
        x = -1.0;   y = -1.0;
        CPPUNIT_ASSERT(vector2NearestNeighbourEuclidean(A, B, C, D, a, b, c, d, x, y) == c);
        x = 1.0;   y = -1.0;
        CPPUNIT_ASSERT(vector2NearestNeighbourEuclidean(A, B, C, D, a, b, c, d, x, y) == d);
    }
    
    /// test bilinearinterpolation
    void testBilinearInterpolation() {
        
        double a = 2.0;
        double b = 2.0;
        double c = 4.0;
        double d = 4.0;
        double x = 0.0;
        double y = 0.0;
        
        double x1 = -2.0;
        double x2 = 2.0;
        double y1 = -2.0;
        double y2 = 2.0;
        
        double val = bilinearInterpolation(x1, x2, y1, y2, a, b, c, d, x, y);
        
        CPPUNIT_ASSERT(abs(val- 3.0) < EPSILON);
        
    }
    
    /// test point in Poly helper function
    void testPointInPoly() {
        Vector2<double> p[4];
        p[3] = double2(1,1);
        p[2] = double2(6,2);
        p[1] = double2(7,5);
        p[0] = double2(0,6);
        
        CPPUNIT_ASSERT(pointInCell(p, double2(3,3)) == true);
        CPPUNIT_ASSERT(pointInCell(p, double2(8,1)) == false);
        CPPUNIT_ASSERT(pointInCell(p, double2(1,5)) == true);
        CPPUNIT_ASSERT(pointInCell(p, double2(3.5,1.5)) == true); // on poly line --> shall return true
        
#warning not yet implemented
        /*
        // new test, inspired of triangle
        p[0] = double2(0.5, 0.5);
        p[1] = double2(1.5, 0.5);
        p[2] = double2(1, 1.5);
        CPPUNIT_ASSERT(pointInPoly(p, 3, double2(1,1)) == true);
        */
        // now test special cases edges
        p[0] = double2(0, 2);
        p[1] = double2(1, 2);
        p[2] = double2(1, 1);
        p[3] = double2(0, 1);
        
        // edge points should be contained!
        for(int i = 0; i < 4; i++) {
          CPPUNIT_ASSERT(pointInCell(p, p[i]));                          // edge points
          CPPUNIT_ASSERT(pointInCell(p, 0.5 * (p[i] + p[(i + 1) % 4])));// middle points on edges
        }
    }
    /// test if haversine formula works correctly
    void testLLDistance() {
        CPPUNIT_ASSERT(fequals(1546.48804834919, distanceLLdeg(double2(10, 10), double2(20, 20), EARTH_RADIUS))); // precalculated with matlab
        CPPUNIT_ASSERT(fequals(0, distanceLLdeg(double2(-10, 30), double2(-10, 30)))); // same point
        
        CPPUNIT_ASSERT(fequals(6581.0845200886, distanceLLdeg(double2(-10, 30), double2(50, 10), EARTH_RADIUS))); // precalculated with matlab
        
        
        // new distance checker!!!
        double2 p1 = double2(10, 10);
        double2 p2 = double2(30, 30);
        
        
        double dist = distanceLLdeg(p1, p2);
        
        double3 p1xyz = spherical2xyz(vector2deg2rad(p1));
        double3 p2xyz = spherical2xyz(vector2deg2rad(p2));
        
        double cosdist = acos(p1xyz * p2xyz);
        std::cout<<"distance haversine: "<<dist<<" distance xy cos"<<cosdist<<std::endl;
        CPPUNIT_ASSERT(fequals(dist, acos(p1xyz * p2xyz)));
        
        
        // sum up distances
        double sum = 0;
        for(int i = 1; i <=360; i++) {
            sum += distanceLLdeg(double2(i - 1, 0), double2(i, 0));
        }
        CPPUNIT_ASSERT(fequals(2 * PI, sum));
    }
    
    
    /// test if triangle area function works properly
    void testLLTriangleAndPoly() {
        double area = triangleLLdegArea(double2(-40, 30), double2(10, 32), double2(-20, 10), EARTH_RADIUS);
        CPPUNIT_ASSERT(fequals(6270020.2190956147388, area));
        
        return;
        
        // test also Poly
        Vector2<double> p[4];
        p[0] = double2(1,1);
        p[1] = double2(6,2);
        p[2] = double2(7,5);
        p[3] = double2(0,6);
        
        CPPUNIT_ASSERT(fequals(297226.119272119598, convexPolygonLLdegArea(p, 4, EARTH_RADIUS))); // precalculated
        
        // test with xy conversion
        Vector3<double> pXYZ[4];
#pragma warning fix here later!
        double offx = 4; // this should be later obsolete!!!
        double offy = 4;
        for(int i = 0; i < 4; i++) {
            pXYZ[i] = spherical2xyz(vector2deg2rad(p[i] + double2(offx, offy)));
        }
        
        CPPUNIT_ASSERT(fequals(297226.119272119598, convexPolygonXYZArea(pXYZ, 4, EARTH_RADIUS))); // precalculated
    }
    
    /// test intersection between quad and triangle
    void testQuadIntersection() {
        // first create ABC and quad
        double2 ABC[3];
        ABC[0] = double2(2, 0);
        ABC[1] = double2(4, 9);
        ABC[2] = double2(6, 0);
        
        double2 Quad[4];
        Quad[0] = double2(2,7);
        Quad[1] = double2(7,8);
        Quad[2] = double2(8,1);
        Quad[3] = double2(1,2);
        
        int count = 0;
        double2 *ipoints = new double2[20];
        CPPUNIT_ASSERT(quadrilateralTriangleIntersection(ABC, Quad, ipoints, &count));
        
        CPPUNIT_ASSERT(count == 4);
        CPPUNIT_ASSERT(valueInArray(ipoints, count, double2(4.34042553191, 7.468085106)));
        CPPUNIT_ASSERT(valueInArray(ipoints, count, double2(3.6279069767441, 7.3255813953488)));
        CPPUNIT_ASSERT(valueInArray(ipoints, count, double2(5.7049180327868, 1.3278688524)));
        CPPUNIT_ASSERT(valueInArray(ipoints, count, double2(2.3999999999, 1.8)));
        
        
        // 2nd test
        ABC[0] = double2(2, 3);
        ABC[1] = double2(5, 7);
        ABC[2] = double2(-0.4, 12);
        
        Quad[0] = double2(0.2,2.5);
        Quad[1] = double2(0.2,3.0);
        Quad[2] = double2(-0.2,3.0);
        Quad[3] = double2(-0.2,2.5);
        
        count = 0;
        CPPUNIT_ASSERT(!quadrilateralTriangleIntersection(ABC, Quad, ipoints, &count));
        CPPUNIT_ASSERT(count == 0);
        
        
        SafeDeleteA(ipoints);
        
    }
    
    /// test if part ratio computations are correct!
    void testPartRatio() {
        double2 A = double2(0, 0);
        double2 B = double2(6, 4);
        double2 T = double2(3, 2);
        
        // should be 1!
        CPPUNIT_ASSERT(fequals(partRatio(A, B, T), 1));
        CPPUNIT_ASSERT(fequals(partRatio(B, A, T), 1));
        CPPUNIT_ASSERT(fequals(partRatio(T, A, B), -0.5));
    }
    
    void testSort3El() {
        // make 3 dummy vectors
        double2 A;
        double2 B;
        double2 C;
        double2 rA = double2(-1, 3);
        double2 rB = double2(0, 5);
        double2 rC = double2(2, 7);
        
        // check all combinations
        A = rA; B = rB; C = rC;                         // ABC
        sortABClon(A, B, C);
        CPPUNIT_ASSERT(A == rA && B == rB && C == rC);
        
        
        A = rA; B = rC; C = rB;                         // ACB
        sortABClon(A, B, C);
        CPPUNIT_ASSERT(A == rA && B == rB && C == rC);
        
        
        A = rB; B = rA; C = rC;                         //BAC
        sortABClon(A, B, C);
        CPPUNIT_ASSERT(A == rA && B == rB && C == rC);
        
        
        A = rB; B = rC; C = rA;                         //BCA
        sortABClon(A, B, C);
        CPPUNIT_ASSERT(A == rA && B == rB && C == rC);
        
        
        A = rC; B = rA; C = rB;                         //CAB
        sortABClon(A, B, C);
        CPPUNIT_ASSERT(A == rA && B == rB && C == rC);
        
        
        A = rC; B = rB; C = rA;                         //CBA
        sortABClon(A, B, C);
        CPPUNIT_ASSERT(A == rA && B == rB && C == rC);
    }
    
    void testIndexTransformations() {
        // we set w = 6, h = 4
        int w = 4;
        int h = 4;
        // start index
        int i = 1;
        int j = 0;
        
        int2 res = -1;

        // a bit more complicated test --> construct therefore an array where each element has a unique number!
        int *a = new int[w * h];
        int *b = new int[w * h];
        for(int k = 0; k < w * h; k++) {
            a[k] = k + 1;
            b[k] = k + 1;
        }
        using namespace std;
        
        // now rotate array 90 CW
        int oldw = w;
        int oldh = h;
        arrayRotate90CW(a, w, h);
        
        cout<<endl<<"original array"<<endl;
        
        printArrayNEMOStyle(b, oldw, oldh);
        cout<<"90° clock wise rotated array"<<endl;
        printArrayNEMOStyle(a, w, h);
        
        // test
        for(int x = 0; x < w; x++)
            for(int y = 0; y < h; y++) {
                res = rotateIndices90CW(x, y, oldw, oldh);
                int valA = a[res.x + res.y * w];
                int valB = b[x + y * oldw];
                CPPUNIT_ASSERT(valA == valB);
            }
        
        // reset
        for(int k = 0; k < w * h; k++) {
            a[k] = k + 1;
            b[k] = k + 1;
        }
        
        // now rotate array 90 CCW
        oldw = w;
        oldh = h;
        arrayRotate90CCW(a, w, h);
        
        cout<<endl<<"original array"<<endl;
        
        printArrayNEMOStyle(b, oldw, oldh);
        cout<<"90° counter clock wise rotated array"<<endl;
        printArrayNEMOStyle(a, w, h);
        
        // test
        for(int x = 0; x < w; x++)
            for(int y = 0; y < h; y++) {
                res = rotateIndices90CCW(x, y, oldw, oldh);
                int valA = a[res.x + res.y * w];
                int valB = b[x + y * oldw];
                CPPUNIT_ASSERT(valA == valB);
            }

        // mirror X coords
        for(int k = 0; k < w * h; k++) {
            a[k] = k + 1;
            b[k] = k + 1;
        }
        
        // now rotate array 90 CCW
        oldw = w;
        oldh = h;
        arrayMirrorX(a, w, h);
        
        cout<<endl<<"original array"<<endl;
        
        printArrayNEMOStyle(b, oldw, oldh);
        cout<<"flipped array"<<endl;
        printArrayNEMOStyle(a, w, h);
        
        // test
        for(int x = 0; x < w; x++)
            for(int y = 0; y < h; y++) {
                res = mirrorIndicesX(x, y, oldw, oldh);
                int valA = a[res.x + res.y * w];
                int valB = b[x + y * oldw];
                CPPUNIT_ASSERT(valA == valB);
            }
        
        SafeDeleteA(a);
        SafeDeleteA(b);
    }
};

#endif
