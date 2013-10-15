//
//  MetricsTest.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 01.08.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef nemoExtractor_MetricsTest_h
#define nemoExtractor_MetricsTest_h

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


#include "Metrics.h"

// see http://cppunit.sourceforge.net/doc/lastest/cppunit_cookbook.html for details

/// Vector2 class tests
class MetricsTest : public CppUnit::TestFixture
{
private:
	//fast setup macros
	CPPUNIT_TEST_SUITE(MetricsTest);
    CPPUNIT_TEST(testOrientation); // orientation handling tests
    CPPUNIT_TEST(testArrayRotation);
    CPPUNIT_TEST(testArrayReorientation);
    // Those two functions here, require the surface integral to be solved more properly!
    //CPPUNIT_TEST(testllTriangle);
    //CPPUNIT_TEST(testllPolygon);
    
    CPPUNIT_TEST(testllVolumeMetricsEasy); // line based tests
    CPPUNIT_TEST(testllAreaMetricsEasy);
    
    CPPUNIT_TEST(testllSplineAreaMetricsEasy); // spline based tests
    CPPUNIT_TEST(testllSplineAreaMetricsCurvy);
    CPPUNIT_TEST(testPolyDerivative);
    CPPUNIT_TEST(testllSplineAreaMetricsConvex);
    CPPUNIT_TEST(testSplineIntegral);
	CPPUNIT_TEST_SUITE_END();
    
    
public:
    
    void testSplineIntegral() {
        double a = 1;
        double b = 6;
        double spline[4] = {4, 3, 2, 1};
        double2 A = double2(0, 0);
        double2 B = double2(5, 0);
        
        double res = Metrics::splineLineIntegral(a, b, spline, A, B);
        double correct = 6475.0 / 12.0;//wolfram alpha
        CPPUNIT_ASSERT(fequals(res, correct));
        
        res = Metrics::splineLineIntegral(b, a, spline, A, B);
        correct = -6475.0 / 12.0;// calculus
        CPPUNIT_ASSERT(fequals(res, correct));
        
    }
    void testPolyDerivative() {
        double poly[6] = {1,3,5,7,-40,1};
        double res[6] = {3,10,21,-160,5, 0};
        Metrics::polyDerivative(poly, 6);
        for(int i = 0; i < 6; i++)CPPUNIT_ASSERT(fequals(poly[i], res[i]));
    }
    
    void testllSplineAreaMetricsConvex() {
        /// test if convex cases are handled propely(no double counting!!!)
        double2 A = double2(2.8, 3.2);
        double2 B = double2(7.3, 3.2);
        
        int w = 4;
        int h = 4;
        double2 *a = new double2[w * h];
        //manual construction
        a[0 + 0 * 4] = double2(3.5, 1); // first row
        a[1 + 0 * 4] = double2(5, 1.1);
        a[2 + 0 * 4] = double2(7, 1.05);
        a[3 + 0 * 4] = double2(7.5, 1);
        
        a[0 + 1 * 4] = double2(2.5, 1.5); // second row
        a[1 + 1 * 4] = double2(4, 1.8);
        a[2 + 1 * 4] = double2(6, 2);
        a[3 + 1 * 4] = double2(8, 1.4);
        
        a[0 + 2 * 4] = double2(1, 2); // third row
        a[1 + 2 * 4] = double2(3, 3);
        a[2 + 2 * 4] = double2(6.5, 3.2);
        a[3 + 2 * 4] = double2(9, 2);
        
        a[0 + 3 * 4] = double2(1.5, 4); // fourth row
        a[1 + 3 * 4] = double2(4.5, 4.1);
        a[2 + 3 * 4] = double2(7.5, 4);
        a[3 + 3 * 4] = double2(10, 4.2);
        
        double *depths = new double[w * h];
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++) {
                depths[i + j * w] = 1.0;
            }
        
        // special case, edge points of cells lie on line
        double correct = 0.001 * Metrics::llDistance(A, B);
        double test = Metrics::splineArea(A, B, a, depths, w, h);
        CPPUNIT_ASSERT(fequals(correct, test));
        
        // special case, vertical line x = 3.2
        A = double2(3.2, 2);
        B = double2(3.2, 3.7);
        correct = 0.001 * Metrics::llDistance(A, B);
        test = Metrics::splineArea(A, B, a, depths, w, h);
        CPPUNIT_ASSERT(fequals(correct, test));
        
        // normal case, random line
        A = double2(2, 2.3);
        B = double2(7, 3.7);
        correct = 0.001 * Metrics::llDistance(A, B);
        test = Metrics::splineArea(A, B, a, depths, w, h);
        CPPUNIT_ASSERT(fequals(correct, test));
        
        // special case, 3 intersections for one line (plot([1.75 4.75], [2.2 4.2], 'g'))
        A = double2(1.75, 2.2);
        B = double2(4.6, 4.05);
        correct = 0.001 * Metrics::llDistance(A, B);
        test = Metrics::splineArea(A, B, a, depths, w, h);
        CPPUNIT_ASSERT(fequals(correct, test));
        
        SafeDeleteA(depths);
        SafeDeleteA(a);
    }
    
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
        
        double2 *a = new double2[w * h];
        double *depths = new double[w * h];
        
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++) {
                // apply some distortion
                double2 v = double2(i * 10, j * 10);
                v = double2(v.x - 0.004 * v.y * v.y, v.y);
                v = double2(v.x, v.y - 0.001 * v.x * v.x + 0.0001 * v.x);
                a[i + j * w] = v;
                
                depths[i + j * w] = 1.0;
            }
        
        double correct = 0.001 * Metrics::llDistance(A, B);
        double test = Metrics::splineArea(A, B, a, depths, w, h);
        
        CPPUNIT_ASSERT(fequals(correct, test));
        
        SafeDeleteA(depths);
        SafeDeleteA(a);
    }
    
    void testArrayRotation() {
        
        using namespace std;
        
        double2 *a = new double2[4];
        int orientation = 0;
        double x = 70.0;
        int w = 2, h = 2;
        // NEMO
        a[0 + 2 * 0] = double2(0, 0);
        a[1 + 2 * 0] = double2(x, 0);
        a[0 + 2 * 1] = double2(0, x);
        a[1 + 2 * 1] = double2(x, x);
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // now rotate clockwise/counterclockwise
        cout<<endl;
        printArrayNEMOStyle(a, w, h);
        arrayRotate90CW(a, w, h);
        orientation = Metrics::orientationOfCoords(a, w, h);
        cout<<endl;
        printArrayNEMOStyle(a, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO90CW);
        
        arrayRotate90CCW(a, w, h);
        orientation = Metrics::orientationOfCoords(a, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        SafeDeleteA(a);
    }
    void testOrientation() {
        using namespace std;
//#error redo test here!
        double2 a[4];
        int orientation = 0;
    
        int w = 2;
        int h = 3;
        double2 *testarray = new double2[w * h];
        double2 *referencearray = new double2[w * h];
        
        // fill reference array with NEMO orientation
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++)
                referencearray[i + j * w] = double2(i, j);
        
        // assert right orientation of reference array
        orientation = Metrics::orientationOfCoords(referencearray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        double x = 70.0;
        
        memcpy(testarray, referencearray, sizeof(double2) * w * h);
        
        
        //NEMO 90CW
        cout<<"before orientation:"<<endl;
        printArrayNEMOStyle(testarray, w, h);
        cout<<"transpose"<<endl;
       // arrayTranspose(testarray, w, h);
        printArrayNEMOStyle(testarray, w, h);
        arrayRotate90CW(testarray, w, h);
        
//#error matrix transpose errorenous!
        cout<<"after orientation:"<<endl;
        printArrayNEMOStyle(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO90CW);
        
        // NEMO90CWMX
        arrayMirrorX(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO90CWMX);
        arrayMirrorX(testarray, w, h); // revert mirroring
        
        
        //NEMO 180CW
        arrayRotate90CW(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO180CW);
        // NEMO180CWMX
        arrayMirrorX(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO180CWMX);
        arrayMirrorX(testarray, w, h); // revert mirroring
        
        
        //NEMO 270CW
        arrayRotate90CW(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO270CW);
        // NEMO270CWMX
        arrayMirrorX(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO270CWMX);
        arrayMirrorX(testarray, w, h); // revert mirroring
        
        //NEMO
        arrayRotate90CW(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        // NEMOMX
        arrayMirrorX(testarray, w, h);
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NEMOMX);
        arrayMirrorX(testarray, w, h); // revert mirroring
        
        
        // now test if invalid data is detected
        testarray[w/2 + h/2 * w] = double2(w * w * w, -h * h * h); // a very big, invalid value!
        orientation = Metrics::orientationOfCoords(testarray, w, h);
        CPPUNIT_ASSERT(orientation == OM_NONE);
    }
    
    void testArrayReorientation() {
        double2 a[4];
        int orientation = 0;
        int w=2, h=2;
        double x = 70.0;
        
        // non oriented data
        a[0 + 2 * 1] = double2(x, 0); a[1 + 2 * 1] = double2(-x, x);
        a[0 + 2 * 0] = double2(0, 0); a[1 + 2 * 0] = double2(x, -x);
        CPPUNIT_ASSERT(!Metrics::reorientArray(a, w, h));
        
        // NEMO
        a[0 + 2 * 0] = double2(0, 0);
        a[1 + 2 * 0] = double2(x, 0);
        a[0 + 2 * 1] = double2(0, x);
        a[1 + 2 * 1] = double2(x, x);
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // NEMO90CW
        a[0 + 2 * 0] = double2(x, 0);
        a[1 + 2 * 0] = double2(x, x);
        a[0 + 2 * 1] = double2(0, 0);
        a[1 + 2 * 1] = double2(0, x);
        
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // NEMO180CW
        a[0 + 2 * 0] = double2(x, x);
        a[1 + 2 * 0] = double2(0, x);
        a[0 + 2 * 1] = double2(x, 0);
        a[1 + 2 * 1] = double2(0, 0);
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // NEMO270CW
        a[0 + 2 * 1] = double2(x, x); a[1 + 2 * 1] = double2(x, 0);
        a[0 + 2 * 0] = double2(0, x); a[1 + 2 * 0] = double2(0, 0);
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // NEMOMX
        a[0 + 2 * 1] = double2(x, x); a[1 + 2 * 1] = double2(0, x);
        a[0 + 2 * 0] = double2(x, 0); a[1 + 2 * 0] = double2(0, 0);
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // NEMO90CWMX
        a[0 + 2 * 1] = double2(0, x); a[1 + 2 * 1] = double2(0, 0);
        a[0 + 2 * 0] = double2(x, x); a[1 + 2 * 0] = double2(x, 0);
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // NEMO180CWMX
        a[0 + 2 * 1] = double2(0, 0); a[1 + 2 * 1] = double2(x, 0);
        a[0 + 2 * 0] = double2(0, x); a[1 + 2 * 0] = double2(x, x);
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
        
        // NEMO270CWMX
        a[0 + 2 * 1] = double2(x, 0); a[1 + 2 * 1] = double2(x, x);
        a[0 + 2 * 0] = double2(0, 0); a[1 + 2 * 0] = double2(0, x);
        CPPUNIT_ASSERT(Metrics::reorientArray(a, w, h));
        orientation = Metrics::orientationOfCoords(a, 2, 2);
        CPPUNIT_ASSERT(orientation == OM_NEMO);
    }
    
    void testllVolumeMetricsEasy() {
        int w = 3;
        int h = 3;
        
        // create some nice testarray
        double2 *points = new double2[w * h];
        double *depths = new double[w * h];
        
        // fill arrays
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++)
            {
                points[i + w * j] = double2(j, i);                depths[i + w * j] = 1;
                
                std::cout<<"("<<i<<","<<j<<") = ["<<points[i + w * j].lon<<","<<points[i + h * j].lat<<"]"<<std::endl;
            }
        double2 A = double2(0.5, 0.5);
        double2 B = double2(1.5, 0.5);
        double2 C = double2(1, 1.5);
        
        // sort points
        sortABClon(A, B, C);
        
        double correctvolume = 0.001 * Metrics::llTriangle(A, B, C);
        double volume = Metrics::linearFrustrum(A, B, C, points, depths, w, h);
        std::cout<<">> difference between exact & computed result: "<<correctvolume-volume<<std::endl;
        CPPUNIT_ASSERT(fequals(volume, correctvolume));
        
        SafeDeleteA(points);
        SafeDeleteA(depths);
    }
    
    void testllAreaMetricsEasy() {
        int w = 3;
        int h = 3;
        
        // create some nice testarray
        double2 *points = new double2[w * h];
        double *depths = new double[w * h];
        
        // fill arrays
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++)
            {
                points[i + w * j] = double2(j, i);
                depths[i + w * j] = 1;
                
                std::cout<<"("<<i<<","<<j<<") = ["<<points[i + w * j].lon<<","<<points[i + h * j].lat<<"]"<<std::endl;
            }
        double2 A = double2(0.5, 0.5);
        double2 B = double2(1.5, 0.5);
        
        double correctarea = 0.001 * Metrics::llDistance(A, B);
        double area = Metrics::linearArea(A, B, points, depths, w, h);
        std::cout<<">> difference between exact & computed result: "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        SafeDeleteA(points);
        SafeDeleteA(depths);
    }
    
    void testllSplineAreaMetricsEasy() {
        int w = 3;
        int h = 3;
        
        double offlon = 0.0;
        double offlat = 0.0;
        //offlon = 15.0;
        //offlat = 24.8;
        offlon = 0;//14.0; // this gives amazing exact results!
        offlat = 0;//24.8;
        
        // create some nice testarray
        double2 *points = new double2[w * h];
        double *depths = new double[w * h];
        
        // fill arrays
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++)
            {
                points[i + w * j] = double2(j +offlon, i+offlat);//vector2deg2rad(double2(j+offlon, i+offlat));
                //points[i + w * j] = double2(-2 +i +offlon,-2+ j+offlat);//vector2deg2rad(double2(j+offlon, i+offlat));
                depths[i + w * j] = 1;
                
                std::cout<<"("<<i<<","<<j<<") = ["<<points[i + w * j].lon<<","<<points[i + h * j].lat<<"]"<<std::endl;
            }
        double2 A = double2(0.5+offlon, 0.5+offlat);//vector2deg2rad(double2(0.5+offlon, 0.5+offlat));
        double2 B = double2(1.5+offlon, 0.5+offlat);//vector2deg2rad(double2(1.5+offlon, 0.5+offlat));
        
        double correctarea = 0.001 * Metrics::llDistance(A, B);
        double area = Metrics::splineArea(A, B, points, depths, w, h);
        std::cout<<">> difference between exact & computed result: "<<correctarea - area<<std::endl;
        CPPUNIT_ASSERT(fequals(area, correctarea));
        
        SafeDeleteA(points);
        SafeDeleteA(depths);
    }
    
    /// new test for the triangle function
    void testllTriangle() {
        
        // this function will fail!!! ignore and keep on moving the implementation!
        /*
         
         // construct 4 edge points
         double2 P1 = double2(-PI, PI/2);
         double2 P2 = double2(PI, PI/2);
         double2 P3 = double2(PI, -PI/2);
         double2 P4 = double2(-PI, -PI/2);
         
         // calc triangle areas
         double A1 = Metrics::llTriangle(P1, P2, P3);
         double A2 = Metrics::llTriangle(P1, P4, P3);
         
         CPPUNIT_ASSERT(fequals(A1 + A2, 4 * PI));
         
         // a bit more complicated test(divide region in 4 axis aligned triangles!)
         double2 PMU = double2(0.4, PI/2);
         double2 PMD = double2(0.4, -PI/2);
         A1 = Metrics::llTriangle(P1, PMU, PMD);
         A2 = Metrics::llTriangle(P1, P4, PMD);
         
         CPPUNIT_ASSERT(fequals(A1 + A2, 2 * PI));
         
         double A3 = Metrics::llTriangle(PMU, PMD, P2);
         double A4 = Metrics::llTriangle(PMU, PMD, P3);
         CPPUNIT_ASSERT(fequals(A3 + A4, 2 * PI));
         
         CPPUNIT_ASSERT(fequals(A1+A2+A3+A4, 4 * PI));
         
         // now take in the middle a non axis aligned triangle and test if sum is still equal to 4pi!
         double2 A = double2(-PI/2, -PI/2);
         double2 B = double2(PI/4, PI/2);
         double2 C = double2(PI, 0);
         
         double A5 = Metrics::llTriangle(A, B, C);
         A2 = Metrics::llTriangle(P1, A, B);
         A1 = Metrics::llTriangle(P1, P4, A);
         A3 = Metrics::llTriangle(B, P2, C);
         A4 = Metrics::llTriangle(A, C, P3);
         
         double A6 = Metrics::llTriangle(P1, B, P2);
         double A7 = Metrics::llTriangle(P4, A, P3);
         
         double sum = A1 + A2 + A3 + A4 + A5 + A6 + A7;
         //CPPUNIT_ASSERT(fequals(sum, 4 * PI));*/
    }
    
    //    void testllPolygon() {
    //        // construct 4 edge points
    //        double2 P1 = double2(-PI, PI/2);
    //        double2 P2 = double2(PI, PI/2);
    //        double2 P3 = double2(PI, -PI/2);
    //        double2 P4 = double2(-PI, -PI/2);
    //
    //        // 4 triangle points
    //        double2 A = double2(-PI/2, PI/4);
    //        double2 B = double2(PI/2, PI/4);
    //        double2 C = double2(0, -PI/4);
    //
    //        double3 p[6];
    //        // now calc data for a triangulation with no axis aligned triangles
    //
    //        p[0] = Hom(P1);
    //        p[1] = Hom(P2);
    //        p[2] = Hom(B);
    //        p[3] = Hom(A);
    //        double A2 = Metrics::llPolygonArea(p, 4);
    //        double test = Metrics::llTriangle(A, B, P2);
    //        double C2 = Metrics::llTriangle(A, B, P2) + Metrics::llTriangle(P1, A, P2);
    //        //CPPUNIT_ASSERT(fequals(A2, C2));
    //
    //        p[0] = Hom(P1);
    //        p[1] = Hom(A);
    //        p[2] = Hom(C);
    //        p[3] = Hom(double2(0, -PI/2));
    //        p[4] = Hom(P4);
    //        double A3 = Metrics::llPolygonArea(p, 5);
    //        double C3 = Metrics::llTriangle(P1, P4, A) + Metrics::llTriangle(P4, A,  C) + Metrics::llTriangle(P4, C, double2(0, -PI/2));
    //       // CPPUNIT_ASSERT(fequals(A3, C3));
    //
    //
    //        p[0] = Hom(P2);
    //        p[1] = Hom(B);
    //        p[2] = Hom(C);
    //        p[3] = Hom(double2(0, -PI/2));
    //        p[4] = Hom(P3);
    //        double C4 = Metrics::llTriangle(double2(0, -PI/2), C, P3) + Metrics::llTriangle(C, B, P3) + Metrics::llTriangle(B, P2, P3);
    //        double A4 = Metrics::llPolygonArea(p, 5);
    //
    //        sortABClon(A, B, C);
    //        // calc triangleArea
    //        double A1 = Metrics::llTriangle(A, B, C);
    //
    //        CPPUNIT_ASSERT(fequals(A1, 2.2214414690791831235079404950303468493073108446878451));
    //
    //        double csum = C2+C3+C4+A1;
    //        double sum = A1+A2+A3+A4;
    //        CPPUNIT_ASSERT(fequals(csum, 4 * PI));
    //        CPPUNIT_ASSERT(fequals(sum, 4 * PI));
    //        
    //    }
};

#endif
