//
//  AreaCalculatorTest.h
//  gmetric
//
//  Created by Leonhard Spiegelberg on 18.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef gmetric_AreaCalculatorTest_h
#define gmetric_AreaCalculatorTest_h

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

#include "AreaCalculator.h"

class AreaCalculatorTest : public CppUnit::TestFixture
{
private:
	//fast setup macros
	CPPUNIT_TEST_SUITE(AreaCalculatorTest);
    CPPUNIT_TEST(testLinear);
    CPPUNIT_TEST(testllSplineAreaMetricsCurvy);
    CPPUNIT_TEST(testllSplineAreaMetricsConvex);
    CPPUNIT_TEST(testDistance);
	CPPUNIT_TEST_SUITE_END();
    
public:
    
    /// test function to test if distances between two points are calculated properly!
    void testDistance() {
        double2 A = double2(-55.55, 53.86);
        double2 B = double2(-48.228, 60.5667);
        
        using namespace std;
        cout<<endl<<"distance AR7W line: "<<AreaCalculator::instance().llDistanceUnitSphere(A, B, true) * EARTH_RADIUS<<"km"<<endl;
        
        A = double2(-180, 0);
        B = double2(180, 0);
        double distance = AreaCalculator::instance().llDistanceUnitSphere(A, B, false);
        cout<<"distance: "<<distance<<endl;;
        CPPUNIT_ASSERT(fequals(distance, 2 * PI));
        
        A = double2(-53.98333, 55.263333);
        B = double2(-54.0156269, 55.232443);
        distance = AreaCalculator::instance().llDistanceUnitSphere(A, B, true);
        cout<<"distance: "<<distance<<endl;
    }
    
    /// test an easy scenario of a small grid with cell eges as lines
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
        
        double2 A, B;
        double area, correctarea;
        bool domain;

        A = double2(0.5, 0.5);
        B = double2(1.5, 0.5);
        
        correctarea  = AreaCalculator::instance().llDistanceUnitSphere(A, B, false) * 0.001;
        area         = AreaCalculator::instance().linear(&grid, A, B, "", false);
        std::cout<<">> difference between exact & computed result: "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        // now test with real data
        // TEST
        // .......
        area = AreaCalculator::instance().linear(&grid, A, B, "", true);
        // note that depth is in m^2! distance should return something in km!
        correctarea = AreaCalculator::instance().llDistanceUnitSphere(A, B, true) * 0.001;
        
        
        std::cout<<">> difference between exact & computed result(real world data): "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        // next test with "axis" aligned line
        A = double2(0, 0);
        B = double2(0, 2);
        area = AreaCalculator::instance().linear(&grid, A, B, "", false);
        correctarea = AreaCalculator::instance().llDistanceUnitSphere(A, B, false) * 0.001;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        // both points lies within one cell
        A = double2(0.2, 0.1);
        B = double2(0.6, 0.5);
        area = AreaCalculator::instance().linear(&grid, A, B, "", false);
        correctarea = AreaCalculator::instance().llDistanceUnitSphere(A, B, false) * 0.001;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        // now check helper function
        A = double2(0, 0);
        B = double2(0, 2); // point inside of domain
        domain = AreaCalculator::instance().pointsInAnyLinearCell(&grid, 1, 1, 3, 3, A, B);
        CPPUNIT_ASSERT(domain);
        A = double2(0, -1); // point outside of domain
        B = double2(0, 3); // point outside of domain
        domain = AreaCalculator::instance().pointsInAnyLinearCell(&grid, 1, 1, 3, 3, A, B);
        CPPUNIT_ASSERT(!domain);
        A = double2(0, -1); // point outside of domain
        B = double2(0, 2);
        domain = AreaCalculator::instance().pointsInAnyLinearCell(&grid, 1, 1, 3, 3, A, B);
        CPPUNIT_ASSERT(!domain);
    }
    
    /// a pretty regular curvy grid, to test the spline algorithm for
    void testllSplineAreaMetricsCurvy() {
        using namespace std;
        //        [x,y] = meshgrid(-80:20:80, -80:20:80);
        //
        //        x = x - 0.004 * y.^2;
        //        y = y - 0.001 * x.^2 + 0.0001 * x;
        
        double2 A = double2(5, 5);
        double2 B = double2(30, 60);
        
        int w = 9;
        int h = 9;
        NCGrid grid;
        
        // init a ncgrid
        grid.init(w, h);
        
        // fill with coords & depths' values
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++) {
                // apply some distortion
                double2 v = double2(i * 10, j * 10);
                v = double2(v.x - 0.004 * v.y * v.y, v.y);
                v = double2(v.x, v.y - 0.001 * v.x * v.x + 0.0001 * v.x);
                grid.setCoordinate(i, j, v);
                grid.setDepth(i, j, 1.0);
            }
        grid.orientation = grid.detectOrientation();
        
        
        double correctarea  = AreaCalculator::instance().llDistanceUnitSphere(A, B, false) * 0.001;
        double area         = AreaCalculator::instance().spline(&grid, A, B, "", false);
        std::cout<<">> difference between exact & computed result: "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        correctarea  = AreaCalculator::instance().llDistanceUnitSphere(A, B, true) * 0.001;
        area         = AreaCalculator::instance().spline(&grid, A, B, "", true);
        std::cout<<">> difference between exact & computed result(real World data): "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));
    }
    
    /// a more complicated test case, for a very irregular curvy grid with a lot of special cases!
    void testllSplineAreaMetricsConvex() {
        /// test if convex cases are handled propely(no double counting!!!)
        double2 A = double2(2.8, 3.2);
        double2 B = double2(7.3, 3.2);
        
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
        
        // special case, edge points of cells lie on line
        double correct = 0.001 * AreaCalculator::instance().llDistanceUnitSphere(A, B, false);
        double test = AreaCalculator::instance().spline(&grid, A, B, "", false);
        CPPUNIT_ASSERT(fequals(correct, test));
        
        // special case, vertical line x = 3.2
        A = double2(3.2, 2);
        B = double2(3.2, 3.7);
        correct = 0.001 * AreaCalculator::instance().llDistanceUnitSphere(A, B, false);
        test = AreaCalculator::instance().spline(&grid, A, B, "", false);
        CPPUNIT_ASSERT(fequals(correct, test));
        
        // normal case, random line
        A = double2(2, 2.3);
        B = double2(7, 3.7);
        correct = 0.001 * AreaCalculator::instance().llDistanceUnitSphere(A, B, false);
        test = AreaCalculator::instance().spline(&grid, A, B, "", false);
        CPPUNIT_ASSERT(fequals(correct, test));
        
        // special case, 3 intersections for one line (plot([1.75 4.75], [2.2 4.2], 'g'))
        A = double2(1.75, 2.2);
        B = double2(4.6, 4.05);
        correct = 0.001 * AreaCalculator::instance().llDistanceUnitSphere(A, B, false);
        test = AreaCalculator::instance().spline(&grid, A, B, "", false);
        CPPUNIT_ASSERT(fequals(correct, test));
    }
    
    
};


#endif
