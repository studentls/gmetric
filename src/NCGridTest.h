//
//  NCGridTest.h
//  gmetric
//
//  Created by Leonhard Spiegelberg on 18.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef gmetric_NCGridTest_h
#define gmetric_NCGridTest_h

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

#include "NCGrid.h"

class NCGridTest : public CppUnit::TestFixture
{
private:
	//fast setup macros
	CPPUNIT_TEST_SUITE(NCGridTest);
    CPPUNIT_TEST(testClassHelpers);
	CPPUNIT_TEST_SUITE_END();
    
public:
    
    
    void testClassHelpers() {
        // define here later maybe some nice test routines for the transformation functions...
    }
};

#endif
