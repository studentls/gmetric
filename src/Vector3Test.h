//
//  Vector3Test.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 31.07.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef nemoExtractor_Vector3Test_h
#define nemoExtractor_Vector3Test_h

#include "Vector3.h"

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




// see http://cppunit.sourceforge.net/doc/lastest/cppunit_cookbook.html for details

/// Vector2 class tests
class Vector3Test : public CppUnit::TestFixture
{
private:
	//fast setup macros
	CPPUNIT_TEST_SUITE(Vector2Test);
	CPPUNIT_TEST(testOperators);
	CPPUNIT_TEST_SUITE_END();
    
public:
    
    /// test +, - , *, /
    void testOperators() {
        Vector3<double> v = Vector3<double>(2.0, 5.0, 1.0);
        Vector3<double> v2 = Vector3<double>(2.0, 5.0, 1.0);
        
        CPPUNIT_ASSERT(v == v2);
        CPPUNIT_ASSERT(v.equals(v2));
        v2.x = 2.0 + EPSILON / 2;
        
        CPPUNIT_ASSERT(v == v2);
        CPPUNIT_ASSERT(v.equals(v2));
        
        CPPUNIT_ASSERT(v + v2 == Vector3<double>(4.0, 10.0, 2.0));
        CPPUNIT_ASSERT(v - v2 == Vector3<double>(0.0, 0.0, 0.0));
        CPPUNIT_ASSERT(5.0 * v == Vector3<double>(10.0, 25.0, 5.0));
        CPPUNIT_ASSERT(v * 5.0 == Vector3<double>(10.0, 25.0, 5.0));
        CPPUNIT_ASSERT(v / 2.0 == Vector3<double>(1.0, 2.5, 0.5));
        
        // length & dot product
        CPPUNIT_ASSERT(abs(v * v - 30.0) < EPSILON);
        v = Vector3<double>(4.0, 4.0, 4.0);
        CPPUNIT_ASSERT(abs(v.length() - sqrt(4.0 * 4.0 + 4.0 * 4.0 + 4.0 * 4.0))<EPSILON);
        
        // cross product
        Vector3<double> a = Vector3<double>(2, 1, -1);
        Vector3<double> b = Vector3<double>(-3, 4, 1);
        CPPUNIT_ASSERT(cross(a, b) == Vector3<double>(5, 1, 11));
        CPPUNIT_ASSERT(cross(b, a) == Vector3<double>(-5, -1, -11));
        CPPUNIT_ASSERT(cross(a, b) == a.cross(b));
        
        CPPUNIT_ASSERT(vector3Min(a,b) == Vector3<double>(-3, 1, -1));
        CPPUNIT_ASSERT(vector3Max(a,b) == Vector3<double>(2, 4, 1));
        
        // divide
        v = double3(12.99, 24.31, 25.2);
        double3 t = v;
        v /= v.z;
        CPPUNIT_ASSERT(v == t / t.z);
    }
};
#endif
