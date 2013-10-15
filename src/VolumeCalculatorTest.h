//
//  VolumeCalculatorTest.h
//  gmetric
//
//  Created by Leonhard Spiegelberg on 04.10.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef gmetric_VolumeCalculatorTest_h
#define gmetric_VolumeCalculatorTest_h

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

#include "VolumeCalculator.h"
#include "Vector2.h"

class VolumeCalculatorTest : public CppUnit::TestFixture
{
private:
	//fast setup macros
	CPPUNIT_TEST_SUITE(VolumeCalculatorTest);
    CPPUNIT_TEST(testPoly);
    CPPUNIT_TEST(testRealWorldTriangle);
    CPPUNIT_TEST(testRealWorldTriangleAndPoly);
    CPPUNIT_TEST(testLinear);
    CPPUNIT_TEST(testTriangleWithSplineSegment);
    CPPUNIT_TEST(testTriangleWithSplineSegmentReal);
    CPPUNIT_TEST(testSpline);
	CPPUNIT_TEST_SUITE_END();
    
public:
    
    /// test if faces can be added linearly
    void testRealWorldTriangle() {
        double2 A = double2(5, 2);
        double2 B = double2(1, 9);
        double2 C = double2(2, -4);
        
        // calc center point
        double2 center = 1.0 / 3.0 * (A + B + C);
        
        VolumeCalculator& vc = VolumeCalculator::instance();
        double correctarea = 0;
        double area = 0;
        
        
        correctarea = vc.llTriangleUnitSphere(A, B, C, true);
        area = 0;
        area += vc.llTriangleUnitSphere(A, B, center, true);
        area += vc.llTriangleUnitSphere(A, C, center, true);
        area += vc.llTriangleUnitSphere(B, C, center, true);
        
        std::cout<<">> difference correct - test: "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        
        // 2nd test with different coords
        A = double2(3.2, 3.7);
        B = double2(3, 3.3);
        C = double2(3.5, 3.7);
        center = 1.0 / 3.0 * (A + B + C);
        correctarea = vc.llTriangleUnitSphere(A, B, C, true);
        
        area = 0;
        area += vc.llTriangleUnitSphere(A, B, center, true);
        area += vc.llTriangleUnitSphere(A, C, center, true);
        area += vc.llTriangleUnitSphere(B, C, center, true);
        
        std::cout<<">> difference correct - test: "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));

    }
    
    /// test if the helper function for the spline works properly with real world integration function!
    void testTriangleWithSplineSegmentReal() {
        
        double spline[] = {0, 0, -1, 1}; // = x^3 - x^2
        
        double2 A = double2(0, 0);
        double2 B = double2(1, 0);
        double2 C = double2(0, -5);
        
        // result should be 1/2 *5*1   + 1/4 - 1/3
        double area = 0;
        double correctarea = 0;
        area += VolumeCalculator::instance().llTriangleWithSplineLon(A, B, C, spline, true);
        correctarea += VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, true);
        // now make a combination of two triangles
        C = double2(0, 5);
        area += VolumeCalculator::instance().llTriangleWithSplineLon(A, B, C, spline, true);
        correctarea += VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, true);
        
        // should be equal to 1 now!
        CPPUNIT_ASSERT(fequals(area, correctarea));
    }

    /// test if the helper function for the spline works properly
    void testTriangleWithSplineSegment() {
        
        double spline[] = {0, 0, rad2deg(-1), rad2deg(1)}; // = x^3 - x^2
        
        double2 A = vector2rad2deg(double2(0, 0));
        double2 B = vector2rad2deg(double2(1, 0));
        double2 C = vector2rad2deg(double2(0, -5));
        
        // result should be 1/2 *5*1   + 1/4 - 1/3
        double area = 0;
        double correctarea = 0;
        area += VolumeCalculator::instance().llTriangleWithSplineLon(A, B, C, spline, false);
        correctarea = 0.5 * 5.0    + 0.25- 1.0/3.0;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        // now make a combination of two triangles
        C = vector2rad2deg(double2(0, 5));
        area += VolumeCalculator::instance().llTriangleWithSplineLon(A, B, C, spline, false);
        
        // should be equal to 1 now!
        CPPUNIT_ASSERT(fequals(area, 5));
    }
    
    /// test if the polygon function works properly
    void testPoly() {
        double2 P1 = double2(-1, 1);
        double2 P2 = double2(1, 1);
        double2 P3 = double2(1, -1);
        double2 P4 = double2(-1, -1);
        
        double2 quad[4];
        quad[0] = P1;
        quad[1] = P2;
        quad[2] = P3;
        quad[3] = P4;
        
        VolumeCalculator& vc = VolumeCalculator::instance();
        double2 center = double2(0,0);
        double area = vc.llPolygonAreaUnitSphere(quad, 4, center, false);
        double correctarea = vc.llTriangleUnitSphere(P1, P2, P3, false) + vc.llTriangleUnitSphere(P1, P3, P4, false);
        
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
    }
    
    /// test if the real world triangle area function retrieves good results
    void testRealWorldTriangleAndPoly() {
        
        // test thing, 4 edge points and then 2 random placed points on the edges
        double2 P1 = double2(-180, 90);
        double2 P2 = double2(180, 90);
        double2 P3 = double2(180, -90);
        double2 P4 = double2(-180, -90);
        
        double2 A = double2(10, 90);
        double2 B = double2(-10, -90);
        
        VolumeCalculator& vc = VolumeCalculator::instance();
        
        double area = 0;
        double correctarea = 0;
        
        // first check if area of P1, P2, P3, P4 equals sphere surface area!
        area = vc.llTriangleUnitSphere(P1, P2, P3, true);
        area += vc.llTriangleUnitSphere(P1, P3, P4, true);
        CPPUNIT_ASSERT(fequals(area, 4 * PI));
        
        area = 0;
        // individual triangle check
        area += vc.llTriangleUnitSphere(P1, P4, B, true);
        area += vc.llTriangleUnitSphere(P1, B, A, true);
        area += vc.llTriangleUnitSphere(B, A, P3, true);
        area += vc.llTriangleUnitSphere(A, P3, P2, true);
        CPPUNIT_ASSERT(fequals(area, 4 * PI));
        
        
        // now with center for other special cases
        area = 0;
        double2 C = vector2deg2rad(double2(10, -10)); // center point!
        area += vc.llTriangleUnitSphere(P1, P2, C, true);
        area += vc.llTriangleUnitSphere(P1, P4, C, true);
        area += vc.llTriangleUnitSphere(P4, P3, C, true);
        area += vc.llTriangleUnitSphere(P3, P2, C, true);
        CPPUNIT_ASSERT(fequals(area, 4 * PI));
        
        // now test if the poly function works correctly
        area = 0;
        double2 quad[4];
        double2 center;
        quad[0] = P1;
        quad[1] = P4;
        quad[2] = A;
        quad[3] = B;
        
        for(int i = 0; i < 4; i++)center += quad[i];
        center /= 4.0;
        
        area += vc.llPolygonAreaUnitSphere(quad, 4, center, true);
        
        CPPUNIT_ASSERT(fequals(area, vc.llTriangleUnitSphere(P1, P4, B, true) + vc.llTriangleUnitSphere(P1, B, A, true)));
        
        quad[0] = P3;
        quad[1] = P2;
        quad[2] = A;
        quad[3] = B;
        for(int i = 0; i < 4; i++)center += quad[i];
        center /= 4.0;
        
        area += vc.llPolygonAreaUnitSphere(quad, 4, center, true);
        CPPUNIT_ASSERT(fequals(area, 4 * PI));
        
    }
    
    /// test the linear volume calculation function
    void testLinear() {
        int w = 3;
        int h = 3;
        NCGrid grid;
        
        // init a ncgrid
        grid.init(w, h);
        
        // fill arrays
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++)
            {
                grid.setCoordinate(i, j, double2(j, i));
                grid.setDepth(i, j, 1);
            }
        // check orientation for the latter algorithm
        grid.orientation = grid.detectOrientation();

        double2 A = double2(0.5, 0.5);
        double2 B = double2(1.5, 0.5);
        double2 C = double2(1, 1.5);
        
        // sort points
        sortABClon(A, B, C);
        
        double correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, false);
        double volume = VolumeCalculator::instance().linear(&grid, A, B, C, "", false);
        std::cout<<">> difference between exact & computed result: "<<correctvolume-volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
        
        // now with realworld data
        correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, true);
        volume = VolumeCalculator::instance().linear(&grid, A, B, C, "", true);
        std::cout<<">> difference between exact & computed result(real world): "<<correctvolume-volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));

    }
    
    /// test spline based volume calculation function
    void testSpline() {
        int w = 4;
        int h = 4;
        NCGrid grid;
        grid.init(w, h);
        
        //manual construction
        grid.setCoordinate(0, 0, double2(3.5, 1)); // first row
        grid.setCoordinate(1, 0, double2(5, 1.1));
        grid.setCoordinate(2, 0, double2(7, 1.05));
        grid.setCoordinate(3, 0, double2(7.5, 1));
        
        grid.setCoordinate(0, 1, double2(2.5, 1.5)); // second row
        grid.setCoordinate(1, 1, double2(4, 1.8));
        grid.setCoordinate(2, 1, double2(6, 2));
        grid.setCoordinate(3, 1, double2(8, 1.4));
        
        grid.setCoordinate(0, 2, double2(1, 2)); // third row
        grid.setCoordinate(1, 2, double2(3, 3));
        grid.setCoordinate(2, 2, double2(6.5, 3.2));
        grid.setCoordinate(3, 2, double2(9, 2));
        
        grid.setCoordinate(0, 3, double2(1.5, 4)); // fourth row
        grid.setCoordinate(1, 3, double2(4.5, 4.1));
        grid.setCoordinate(2, 3, double2(7.5, 4));
        grid.setCoordinate(3, 3, double2(10, 4.2));
        
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++) {
                grid.setDepth(i, j, 1.0);
            }
        
        // set orientation
        grid.orientation = grid.detectOrientation();
        
        double2 A, B, C;
        double correctvolume;
        double volume;
        
        // first test case: only 1 triangle point in a cell
        // now the triangle
        A = double2(3.2, 3.7);
        B = double2(2, 2.3);
        C = double2(7, 3.7);
        
        correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, false);
        volume = VolumeCalculator::instance().spline(&grid, A, B, C, "", false);
        std::cout<<">> difference between exact & computed result(spline): "<<correctvolume - volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
        
        
        correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, true);
        volume = VolumeCalculator::instance().spline(&grid, A, B, C, "", true);
        std::cout<<">> difference between exact & computed result(spline, real world): "<<correctvolume - volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
        
        // next testcase: 2 triangle points in one cell!
        A = double2(3.2, 3.7);
        B = double2(1.7, 2.3);
        C = double2(7, 3.7);

        correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, false);
        volume = VolumeCalculator::instance().spline(&grid, A, B, C, "", false);
        std::cout<<">> difference between exact & computed result(spline): "<<correctvolume - volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
        
        correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, true);
        volume = VolumeCalculator::instance().spline(&grid, A, B, C, "", true);
        std::cout<<">> difference between exact & computed result(spline): "<<correctvolume - volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
        
        // next testcase: 3 triangle points in one cell!
        A = double2(3.2, 3.7);
        B = double2(3, 3.3);
        C = double2(3.5, 3.7);
        
        correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, false);
        volume = VolumeCalculator::instance().spline(&grid, A, B, C, "", false);
        std::cout<<">> difference between exact & computed result(spline): "<<correctvolume - volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
        
        correctvolume = 0.001 * VolumeCalculator::instance().llTriangleUnitSphere(A, B, C, true);
        volume = VolumeCalculator::instance().spline(&grid, A, B, C, "", true);
        std::cout<<">> difference between exact & computed result(spline): "<<correctvolume - volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
    }
};

#endif
