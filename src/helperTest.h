//
//  helperTest.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 24.05.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef nemoExtractor_helperTest_h
#define nemoExtractor_helperTest_h

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

/// test of helpers
class helperTest : public CppUnit::TestFixture
{
private:
	//fast setup macros
	CPPUNIT_TEST_SUITE(helperTest);
	CPPUNIT_TEST(testSmallFunctions);
	CPPUNIT_TEST(testQuickSelect);
    CPPUNIT_TEST(testSphericalConversion);
    CPPUNIT_TEST(testCubicSpline);
    CPPUNIT_TEST(testHorner);
    CPPUNIT_TEST(testAitkenNeville);
    CPPUNIT_TEST(testQuadraticRoots);
    CPPUNIT_TEST(testCubicRoots);
	CPPUNIT_TEST_SUITE_END();
public:
    
    void testSmallFunctions() {
        float a = 4.0;
        float b = 6.0;
        CPPUNIT_ASSERT(::min(a,b) == a);
        CPPUNIT_ASSERT(::max(a,b) == b);
    }
    
    void testHorner() {
        double a[5] = {10, 3, -5, 2, 3.4};
        
        double res = evalPolynomial(-2.0, a, 5);
        CPPUNIT_ASSERT(fequals(res, 22.4));
    }
    
    void testQuickSelect() {
        int array[] = {200, 300, 100, 20, 60, 90};
        
	std::cout<<quick_select(array, 0, 6, 3);
        //CPPUNIT_ASSERT(quick_select(array, 0, 6, 3) == 90);
        //CPPUNIT_ASSERT(quick_select(array, 0, 6, 6) == 300);
    }
    
    void testQuadraticRoots() {
        
        cdouble roots[2];
        int numRoots = 0;
        
        // x^2 + 1 = (x-i)(x+i)
        squareRootsPQ(0.0, 1.0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 2);
        CPPUNIT_ASSERT(roots[0] == cdouble(0, 1) || roots[0] == cdouble(0, -1));
        CPPUNIT_ASSERT(roots[1] == cdouble(0, 1) || roots[1] == cdouble(0, -1));
        
        // x^2 - 23x + 112 = (x-16)(x-7)
        squareRootsPQ(-23.0, 112.0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 2);
        CPPUNIT_ASSERT(roots[0] == cdouble(16, 0) || roots[0] == cdouble(7, 0));
        CPPUNIT_ASSERT(roots[1] == cdouble(16, 0) || roots[1] == cdouble(7, 0));
        
        // x^2 +2x+1 = (x+1)(x+1)
        squareRootsPQ(2.0, 1.0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 1);
        CPPUNIT_ASSERT(roots[0] == cdouble(-1, 0));
        
    }
    
    void testCubicRoots() {
        double roots[3];
        int numRoots = 0;
        
        // 3 different real roots, (x-9)(x+3)(x-2) = 54 - 15x -8x^2 + x^3
        cubicRealRoots(1.0, -8.0, -15.0, 54.0, roots, &numRoots);
        
        CPPUNIT_ASSERT(numRoots == 3);
        for(int i = 0; i < 3; i++)CPPUNIT_ASSERT(fequals(roots[i], 9.0) || fequals(roots[i], -3.0) || fequals(roots[i], 2.0));
        
        // 1 real root, -189+27 x-21 x^2+3 x^3 = 3(x-7)(x^2+9)
        cubicRealRoots(3.0, -21.0, 27.0, -189.0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 1);
        CPPUNIT_ASSERT(fequals(7.0, roots[0]));
        
        // special cases, quadratic polynom
        // x^2 +2x+1 = (x+1)(x+1)
        cubicRealRoots(0.0, 1.0, 2.0, 1.0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 1);
        CPPUNIT_ASSERT(fequals(-1.0, roots[0]));
        
        // x^2 - 23x + 112 = (x-16)(x-7)
        cubicRealRoots(0.0, 1.0, -23.0, 112.0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 2);
        CPPUNIT_ASSERT(fequals(roots[0], 16.0) || fequals(roots[0], 7.0));
        CPPUNIT_ASSERT(fequals(roots[1], 16.0) || fequals(roots[1], 7.0));
        
        // x^2 + 1 = (x-i)(x+i)
        cubicRealRoots(0.0, 1.0, 0, 1.0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 0);
        
        // special case linear function
        // 2x+1
        cubicRealRoots(0, 0, 2, 1, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 1);
        CPPUNIT_ASSERT(fequals(roots[0], -0.5));
                       
        // constant function
        // x = 5
        cubicRealRoots(0, 0, 0, 5, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 0);
        // x = 0
        cubicRealRoots(0, 0, 0, 0, roots, &numRoots);
        CPPUNIT_ASSERT(numRoots == 1);
        CPPUNIT_ASSERT(fequals(roots[0], 0.0));
    }
    
    /// test if coordinates are converted correctly
    void testSphericalConversion() {
        double2 p = vector2deg2rad(double2(25, -30));
        
        double3 v = spherical2xyz(p);
        
        CPPUNIT_ASSERT(xyz2spherical(v) == p);
    }
    
    /// test cubic spline
    void testCubicSpline() {
        double x0 = 0;
        double x1 = 1;
        double x2 = 2;
        double x3 = 3;
        double z0 = 4;
        double z1 = 2;
        double z2 = 1;
        double z3 = 3;
        
        double x = 3.0/2.0;
        
        double r = cubicSpline(x, x0, x1, x2, x3, z0, z1, z2, z3);
        CPPUNIT_ASSERT(fequals(r, 5.0 / 4.0));
        
        r = cubicSpline(x1, x0, x1, x2, x3, z0, z1, z2, z3);
        CPPUNIT_ASSERT(fequals(r, z1));
        
        r = cubicSpline(x2, x0, x1, x2, x3, z0, z1, z2, z3);
        CPPUNIT_ASSERT(fequals(r, z2));
        
        
        // now test if coefficients are retrieved correctly
        double a[4];
        cubicSplineCoefficients(x0, x1, x2, x3, z0, z1, z2, z3, a);
        r = evalPolynomial(x1, a, 4);
        CPPUNIT_ASSERT(fequals(r, z1));
        r = evalPolynomial(x2, a, 4);
        CPPUNIT_ASSERT(fequals(r, z2));
        r = evalPolynomial(x, a, 4);
        CPPUNIT_ASSERT(fequals(r, 5.0 / 4.0));
        
        // second test for coefficients, what happens if all x
    }
    
    /// test Aitken Neville Algorithm!
    void testAitkenNeville() {
        double x[4] = {-2, -1, 0, 2};
        double y[4] = {-7, 0, 1, 9};
        double xq = 1;
        
        double res = AitkenNeville(xq, x, y, 4);
        CPPUNIT_ASSERT(fequals(res, 2));
        
        res = AitkenNevilleInplace(xq, x, y, 4);
        CPPUNIT_ASSERT(fequals(res, 2));
    }
};


#endif
