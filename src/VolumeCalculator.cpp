//
//  VolumeCalculator.cpp
//  gmetric
//
//  Created by Leonhard Spiegelberg on 30.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#include "VolumeCalculator.h"

#define BUFFER_SIZE 40

double VolumeCalculator::linear(NCGrid *grid, const double2 _A, const double2 _B, const double2 _C, const std::string filename, const bool realWorldData, const double R) {
    using namespace std;
    
    // (1) first step is to get bounding box of triangle in lon/lat coords
    double2 minABCll = vector2Min(vector2Min(_A, _B), _C);
    double2 maxABCll = vector2Max(vector2Max(_A, _B), _C);
    
    // (2) find bounding box indices!
    int si=-1,ei=-1,sj=-1,ej=-1;
    AreaCalculator::instance().findBoundingBoxIndices(grid, minABCll, maxABCll, si, ei, sj, ej);
    
    // (3) now construct the points boundingBox
    //
    double2 pointsMin = grid->getCoordinate(si, sj);//_points[si + sj * w];
    double2 pointsMax = grid->getCoordinate(ei, ej);//_points[ei + ej * w];
    
    // assert this is a real bounding Box!
    assert(pointsMin.lon <= pointsMax.lon);
    assert(pointsMin.lat <= pointsMax.lat);

    // add triangle points to array
    double2 trianglePoints[] = {_A, _B, _C};
    
    // the volume...
    double volume = 0;
    
    double2 p[4]; // edge points in XYZ
    
    double2 *areaPoints = new double2[BUFFER_SIZE];
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
            p[0] = grid->getCoordinate(i - 1, j);     //points[(i - 1)  + j * w]; // top left
            p[1] = grid->getCoordinate(i,     j);     //points[i        + j * w]; // top right
            p[2] = grid->getCoordinate(i,     j - 1); //points[i        + (j - 1) * w]; // bottom right
            p[3] = grid->getCoordinate(i - 1, j - 1);   //points[(i - 1)  + (j - 1) * w]; // bottom left
            
            // set back to zero
            areaPointsCount = 0;
            
            // now find out which points of A, B, C lie in Cell and vice versa
            for(int k = 0; k < 3; k++) {
                if(pointInCell(p, trianglePoints[k]))areaPoints[areaPointsCount++] = trianglePoints[k];
            }
            for(int k = 0; k < 4; k++) {
                if(pointInTriangle(trianglePoints, p[k]))areaPoints[areaPointsCount++] = p[k];
            }
            
            // now the most interesting part! determine intersection points and add to array
            quadrilateralTriangleIntersection(trianglePoints, p, areaPoints + areaPointsCount, &areaPointsCount);
            
            // if less than 3 points skip! (would mean points are on one line!)
            if(areaPointsCount < 3)continue;
            
            // center point is center of all edge points
            double2 center;
            for(int k = 0; k < areaPointsCount; k++)center += areaPoints[k];
            center *= 1.0 / areaPointsCount;
            
            // get area of spherical polygon, defined by points in vector areaPoints!
            double area = llPolygonAreaUnitSphere(areaPoints, areaPointsCount, center, realWorldData);
            
            assert(area >= 0);
            volume += area * grid->getDepth(i, j);
        }
    
    SafeDeleteA(areaPoints);
    
    // multiply with 0.001 for m to km conversion and R^2 to scale volume from unit sphere to arbitrary sphere with radius R!
    return volume * R * R * 0.001;
}

/// helper function to calc area of a triangle on sphere
/// important! accepts degrees!
double VolumeCalculator::llTriangleUnitSphere(const double2 &_A, const double2 &_B, const double2 &_C, const bool realWorldData) {
    
    // convert values here to radian!
    double2 A = vector2deg2rad(_A);
    double2 B = vector2deg2rad(_B);
    double2 C = vector2deg2rad(_C);
    
    if(!realWorldData) {
        // easy, use cross product for area
        double2 AB = B - A;
        double2 AC = C - A;
        double3 crossproduct = cross(double3(AB.x, AB.y, 0), double3(AC.x, AC.y, 0));
        
        return 0.5 * vector3Length(crossproduct);
    }
    else {
        // now integration over unit sphere again
        
        //(1) first sort after lon
        
        // all three equal?
        if(A == B && B == C)return 0;
        
        sortABClon(A, B, C);
        
        // assert check
        assert(A.lon <= B.lon);
        assert(B.lon <= C.lon);
        
        // handle special cases
        if(fequals(A.lon, B.lon)) {
            // A & B are axis aligned vertically
            
            // easier function
            
            // swap such that A.lat < B.lat!
            if(B.lat < A.lat)swap(A, B);
            
            
            // triangle looks alike
            //    B\
            //    |  \
            //    A----C
            
            double l2 = C.lon;
            double l1 = A.lon;
            
            double b1 = (C.lat - A.lat) / (C.lon - A.lon);
            double b0 = A.lat - b1 * A.lon;
            
            double a1 = (C.lat - B.lat) / (C.lon - B.lon);
            double a0 = B.lat - a1 * B.lon;
            
            // cases...
            double area = 0;
            if(fequals(b1, 0)) {
                area += sin(b0) * (l2 - l1);
            }else {
                area += 1/b1 * (cos(b1 * l2 + b0) - cos(b1 * l1 + b0));
            }
            if(fequals(a1, 0)) {
                area += sin(a0) * (l2 - l1);
            }
            else {
                area += 1/a1 * (cos(a1 * l1 + a0) - cos(a1 * l2 + a0));
            }
            
            return abs(area);
        }
        else if(fequals(B.lon, C.lon)) {
            // B & C are axis aligned vertically
            
            // swap such that B.lat > C.lat
            if(B.lat < C.lat)swap(B, C);
            
            // triangle looks alike
            //        /B
            //      /  |
            //    A----C
            
            double l2 = C.lon;
            double l1 = A.lon;
            
            double b1 = (C.lat - A.lat) / (C.lon - A.lon);
            double b0 = A.lat - b1 * A.lon;
            
            double a1 = (B.lat - A.lat) / (B.lon - A.lon);
            double a0 = A.lat - a1 * A.lon;
            
            // cases...
            double area = 0;
            if(fequals(b1, 0)) {
                area += sin(b0) * (l2 - l1);
            }else {
                area += 1/b1 * (cos(b1 * l2 + b0) - cos(b1 * l1 + b0));
            }
            if(fequals(a1, 0)) {
                area += sin(a0) * (l2 - l1);
            }
            else {
                area += 1/a1 * (cos(a1 * l1 + a0) - cos(a1 * l2 + a0));
            }
            
            return abs(area);
        }
        else {
        
            // no sides are vertically axis aligned
            // meaning division can be performed without harm!
            double b1 = (C.lat - A.lat) / (C.lon - A.lon);
            double b0 = A.lat - b1 * A.lon;
            
            double a01 = (B.lat - A.lat) / (B.lon - A.lon);
            double a00 = A.lat - a01 * A.lon;
            
            double a11 = (C.lat - B.lat) / (C.lon - B.lon);
            double a10 = B.lat - a11 * B.lon;
        
            double area = 0;
            
            double l1 = A.lon;
            double l2 = B.lon;
            double l3 = C.lon;
            
            // cases...
            if(fequals(b1, 0)) {
                area += sin(b0) *(l3 - l1);
            }
            else {
                area += 1 / b1 * (cos(b1 * l1 + b0) - cos(b1 * l3 + b0));
            }
            
            if(fequals(a01, 0)) {
                area += sin(a00) * (l1 - l2);
            }
            else {
                area += 1 / a01 * (cos(a01 * l2 + a00) - cos(a01 * l1 + a00));
            }
            
            if(fequals(a11, 0)) {
                area += sin(a10) * (l2 - l3);
            }
            else {
                area += 1 / a11 * (cos(a11 * l3 + a10) - cos(a11 * l2 + a10));
            }
            
            //
            return abs(area);
        }
    }
}


/// helper function which uses lltriangle to decompose a (convex!) polygon into several small triangles
double VolumeCalculator::llPolygonAreaUnitSphere(double2 *points, const int count, const double2& center, const bool realWorldData) {
    
#warning make here an assert check maybe for the center point! (must lie in polygon!)
    
    // make easy triangulation
    // no area for less than 3 points
    if(count < 3)return 0.0;
    
    // quicker round for 3 points
    if(count == 3)
    {
        double2 tA = points[0];
        double2 tB = points[1];
        double2 tC = points[2];
        return llTriangleUnitSphere(tA, tB, tC, realWorldData);
    }
    
    // first step is to sort points
    // therefore go through points and calc their angle in relation to y axes! (because it can be easily calculated with tan!)

//    // calculate angles
//    double *angles = new double[count];
//    double3 m = cross(center, points[0]);
//    double3 M = cross(double3(0, 0, 1), m); // inf point!
//    M = vector3Normalize(M);
//    for(int i = 1; i < count; i++) {
//        
//        double3 l = cross(center, points[i]);
//        
//        // to retrieve angle, get inf points!
//        double3 L = cross(double3(0, 0, 1), l);
//        
//        // normalize both vectors
//        L = vector3Normalize(L);
//        
//        // set angle, based on dot product of normals!
//        angles[i] = (L * M);
//        
//        assert(-PI <= angles[i] && angles[i] <= PI);
//    }
    
    double *angles = new double[count]; // holds angles between center-point and first point-center
    // calc angles using atan2 to get angles for 4 quadrants!
    for(int i = 0; i < count; i++) {
        double2 L = points[i] - center;
        
        // atan2(y, x) !
        angles[i] = atan2(L.y, L.x);
        
    }
    
    // sort with insertion sort
    for(int i = 1; i < count; i++) {
        for(int k = i; k > 0 && angles[k] < angles[k - 1]; k--) {
            swap(points[k], points[k - 1]);
            swap(angles[k], angles[k - 1]);
        }
    }
    
    int pointcount = count;
    
    double area = 0;
    double2 A = center;
    double2 B;
    double2 C;
    // now calc area based on triangulation
    // A remains fixed as 
    for(int i = 0; i < pointcount; i++) {
        B = points[i];
        C = points[(i + 1) % pointcount];
        
        area += llTriangleUnitSphere(A, B, C, realWorldData);//Metrics::llTriangle(tA, tB, tC, R);
    }
    
    SafeDeleteA(angles);
    
    return area;
}


double VolumeCalculator::llTriangleWithSplineLon(const double2 &splineA, const double2 &splineB, const double2 &C, const double *spline, const bool realWorldData) {
    
    // assert valid pointer
    assert(spline);
    
    // first assert that splineA and splineB have different lon values
    assert(!fequals(splineA.lon, splineB.lon));
    
    // assert order of points!
    assert(splineA.x < splineB.x);
    
    double area = 0;
    
    
    // first add to area, the area of the triangle encaged by splineA, splineB & splineC!
    area += llTriangleUnitSphere(splineA, splineB, C, realWorldData);
    
    
    // simple integration!
    
    // use function of area calculator!
    // for correct results, it is important to check where point C lies!
    // if C lies below the spline segment add directly,
    // iff not change sign!
    double a1 = (splineA.lat - splineB.lat) / (splineA.lon - splineB.lon);
    double a0 = splineA.lat - a1 * splineA.lon;
    if(C.lat < C.lon * a1 + a0) {
        
        area += splineLineIntegral(splineA.x, splineB.x, spline, splineA, splineB, realWorldData);
    }
    else {
        area += -splineLineIntegral(splineA.x, splineB.x, spline, splineA, splineB, realWorldData);
    }
    
    
//#warning test needed here if there is any difference if the order of splineA and splineB is changed!
    return area;
}

/// lat function is the same like lon function, just need to swap coords!
double VolumeCalculator::llTriangleWithSplineLat(const double2 &splineA, const double2 &splineB, const double2 &C, const double *spline, const bool realWorldData) {
    return -llTriangleWithSplineLon(swapXY(splineA), swapXY(splineB), swapXY(C), spline, realWorldData);
}


/// calculate volume using spline based cells
double VolumeCalculator::spline(NCGrid *grid, const double2 _A, const double2 _B, const double2 _C, const std::string filename, const bool realWorldData, const double R) {
    
    using namespace std;
    
    // (1) first step is to get bounding box of triangle in lon/lat coords
    double2 minABCll = vector2Min(vector2Min(_A, _B), _C);
    double2 maxABCll = vector2Max(vector2Max(_A, _B), _C);
    
    // (2) find bounding box indices!
    int si=-1,ei=-1,sj=-1,ej=-1;
    AreaCalculator::instance().findBoundingBoxIndices(grid, minABCll, maxABCll, si, ei, sj, ej);
    
    // (3) now construct the points boundingBox
    //
    double2 pointsMin = grid->getCoordinate(si, sj);//_points[si + sj * w];
    double2 pointsMax = grid->getCoordinate(ei, ej);//_points[ei + ej * w];
    
    // assert this is a real bounding Box!
    assert(pointsMin.lon <= pointsMax.lon);
    assert(pointsMin.lat <= pointsMax.lat);
    
    // add triangle points to array
    double2 trianglePoints[] = {_A, _B, _C};
    
#ifdef MATLAB
    ofs<<"% generated code"<<endl;
    ofs.precision(4);
    ofs<<"fig = figure()"<<endl;
    ofs<<endl;
    for(int i = 0; i < 3; i++) {
        double3 col = double3(0.7, 0.05, 0.05);
        ofs<<"scatter("<<trianglePoints[i].x<<","<<trianglePoints[i].y<<", 15, 'o', 'fill', 'MarkerFaceColor',["<<col.x<<" "<<col.y<<" "<<col.z<<"], 'MarkerEdgeColor', 'none'); hold on;\n";
        ofs<<"plot(["<<trianglePoints[i].x<<" "<<trianglePoints[(i+1) % 3].x<<"],["<<trianglePoints[i].y<<" "<<trianglePoints[(i+1) % 3].y<<"], 'Color',["<<col.x<<" "<<col.y<<" "<<col.z<<"]); hold on;"<<endl;
    }

    ofs<<endl;
#endif
    
    // the volume...
    double volume = 0;
    
    double2 p[4]; // edge points in XYZ
    
    SIntersectionPoint<double> *areaPoints = new SIntersectionPoint<double>[BUFFER_SIZE];
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
            p[0] = grid->getCoordinate(i - 1, j);     //points[(i - 1)  + j * w]; // top left
            p[1] = grid->getCoordinate(i,     j);     //points[i        + j * w]; // top right
            p[2] = grid->getCoordinate(i,     j - 1); //points[i        + (j - 1) * w]; // bottom right
            p[3] = grid->getCoordinate(i - 1, j - 1);   //points[(i - 1)  + (j - 1) * w]; // bottom left
            
            
            
            // construct the four splines describing the cell
            double yLeftSpline[4]; // function lat -> lon
            double yRightSpline[4]; // function lat -> lon
            double xTopSpline[4]; // function lon->lat
            double xBottomSpline[4]; // function lon->lat
            double2 s0, s1, s2, s3; // dummys
            
            /// normal orientation(attention! order is important!!!)
            // left side
            s0 = grid->getCoordinate(i - 1, j - 2);//points.get(i - 1, j - 2);
            s1 = grid->getCoordinate(i - 1, j - 1);//points.get(i - 1, j - 1);
            s2 = grid->getCoordinate(i - 1, j + 0);//points.get(i - 1, j + 0);
            s3 = grid->getCoordinate(i - 1, j + 1);//points.get(i - 1, j + 1);
            cubicSplineCoefficients(s0.y, s1.y, s2.y, s3.y, s0.x, s1.x, s2.x, s3.x, yLeftSpline);
            
            // right side
            s0 = grid->getCoordinate(i + 0, j - 2);//points.get(i + 0, j - 2);
            s1 = grid->getCoordinate(i + 0, j - 1);//points.get(i + 0, j - 1);
            s2 = grid->getCoordinate(i + 0, j + 0);//points.get(i + 0, j + 0);
            s3 = grid->getCoordinate(i + 0, j + 1);//points.get(i + 0, j + 1);
            cubicSplineCoefficients(s0.y, s1.y, s2.y, s3.y, s0.x, s1.x, s2.x, s3.x, yRightSpline);
            
            // top side
            s0 = grid->getCoordinate(i - 2, j);//points.get(i - 2, j);
            s1 = grid->getCoordinate(i - 1, j);//points.get(i - 1, j);
            s2 = grid->getCoordinate(i - 0, j);//points.get(i - 0, j);
            s3 = grid->getCoordinate(i + 1, j);//points.get(i + 1, j);
            cubicSplineCoefficients(s0.x, s1.x, s2.x, s3.x, s0.y, s1.y, s2.y, s3.y, xTopSpline);
            
            // bottom side
            s0 = grid->getCoordinate(i - 2, j - 1);//points.get(i - 2, j - 1);
            s1 = grid->getCoordinate(i - 1, j - 1);//points.get(i - 1, j - 1);
            s2 = grid->getCoordinate(i - 0, j - 1);//points.get(i - 0, j - 1);
            s3 = grid->getCoordinate(i + 1, j - 1);//points.get(i + 1, j - 1);
            cubicSplineCoefficients(s0.x, s1.x, s2.x, s3.x, s0.y, s1.y, s2.y, s3.y, xBottomSpline);

#ifdef MATLAB
            // print spline codes to matlab file
            ofs<<AreaCalculator::instance().printOutSplineCode(xTopSpline, p[0].x, p[1].x, 'x');
            ofs<<AreaCalculator::instance().printOutSplineCode(xBottomSpline, p[3].x, p[2].x, 'x');
            ofs<<AreaCalculator::instance().printOutSplineCode(yLeftSpline, p[0].y, p[3].y, 'y');
            ofs<<AreaCalculator::instance().printOutSplineCode(yRightSpline, p[1].y, p[2].y, 'y');
#endif
            
            // set back to zero
            areaPointsCount = 0;
            
            // now find out which points of A, B, C lie in Cell and vice versa
            for(int k = 0; k < 3; k++) {
               //if(pointInCell(p, trianglePoints[k])) {
                if(AreaCalculator::instance().llPointInSplineCell(trianglePoints[k], yLeftSpline, yRightSpline, xTopSpline, xBottomSpline)) {
                    SIntersectionPoint<double> t;
                    t.point = trianglePoints[k];
                    t.side = SIDE_TRIANGLE;
                    areaPoints[areaPointsCount++] = t;
                }
            }
            
            // check now for all edge points, if they lie within the triangle
            // if so add two points to array with the according sides
            for(int k = 0; k < 4; k++) {
                if(pointInTriangle(trianglePoints, p[k])) {
                    SIntersectionPoint<double> t1, t2;
                    t1.point = p[k];
                    t2.point = p[k];
                    
                    // determine side based on k
                    if(k == 0) {t1.side = SIDE_TOP; t2.side = SIDE_LEFT;} // top left
                    if(k == 1) {t1.side = SIDE_TOP; t2.side = SIDE_RIGHT;} // top right
                    if(k == 2) {t1.side = SIDE_BOTTOM; t2.side = SIDE_RIGHT;} // bottom right
                    if(k == 3) {t1.side = SIDE_BOTTOM; t2.side = SIDE_LEFT;} // bottom left
                    
                    // add both points to array!
                    areaPoints[areaPointsCount++] = t1;
                    areaPoints[areaPointsCount++] = t2;
                }
            }
            
            
            // next step is to intersect all sides with all lines of the triangle and add intersection segments to array!
            intersectTriangleSidesWithCellSides(_A, _B, _C, p, xTopSpline, xBottomSpline, yLeftSpline, yRightSpline, areaPoints, &areaPointsCount);
            
            // if less than 3 points skip! (would mean points are on one line!)
            if(areaPointsCount < 3)continue;
            
            // center point is center of all edge points
            double2 center;
            for(int k = 0; k < areaPointsCount; k++)center += areaPoints[k].point;
            center *= 1.0 / areaPointsCount;
            
            // get area of spherical polygon, defined by points in vector areaPoints!
            double area = llPolygonWithSplinesAreaUnitSphere(areaPoints, areaPointsCount, center, xTopSpline, xBottomSpline, yLeftSpline, yRightSpline, realWorldData);
            
            assert(area >= 0);
            volume += area * grid->getDepth(i, j);
        }
    
    SafeDeleteA(areaPoints);
    
    // multiply with 0.001 for m to km conversion and R^2 to scale volume from unit sphere to arbitrary sphere with radius R!
    return volume * R * R * 0.001;

}

/// function to get area points through intersection routine
void VolumeCalculator::intersectTriangleSidesWithCellSides(const double2 &A, const double2 &B, const double2 &C, const double2 *p, const double *topSpline, const double *bottomSpline, const double *leftSpline, const double *rightSpline, SIntersectionPoint<double> *out, int *count) {
    
    // this function is basically a wrap up to pass different sides to the specified function
    intersectLineSegmentWithCellSides(A, B, p, topSpline, bottomSpline, leftSpline, rightSpline, out, count);
    intersectLineSegmentWithCellSides(A, C, p, topSpline, bottomSpline, leftSpline, rightSpline, out, count);
    intersectLineSegmentWithCellSides(B, C, p, topSpline, bottomSpline, leftSpline, rightSpline, out, count);
}

/// here happens all the magic!
void VolumeCalculator::intersectLineSegmentWithCellSides(const double2 &_A, const double2 &_B, const double2 *p, const double *topSpline, const double *bottomSpline, const double *leftSpline, const double *rightSpline, SIntersectionPoint<double> *out, int *count) {
    // first assert all pointers
    assert(p);
    assert(topSpline);
    assert(bottomSpline);
    assert(leftSpline);
    assert(rightSpline);
    assert(out);
    assert(count);
    
    // next swap A,B such that A.lon < B.lon!
    double2 A = _A;
    double2 B = _B;
    if(A.lon > B.lon)swap(A, B);
    
    assert(A.lon <= B.lon);
    
    // bad hack, but works!
    double2 topleft     = p[0];
    double2 topright    = p[1];
    double2 bottomright = p[2];
    double2 bottomleft  = p[3];
    
#warning actually p0-p4 should be asserted...
    
    // intersect with splines
    // all in all
    double3             intersectionPoints[4 * 3]; // 0-2 top, 3-5 bottom, 6-8 left, 9-11 right
    unsigned int        intersectionPointsCount[4] = {0, 0, 0, 0}; // 0 top, 1 bottom, 2 left, 3 right
    
    AreaCalculator& ac =  AreaCalculator::instance();
    ac.splineLineIntersectionLon(A, B, topSpline, topleft.x, topright.x,
                              intersectionPoints + SIDE_TOP * 3, &intersectionPointsCount[SIDE_TOP]);
    ac.splineLineIntersectionLon(A, B, bottomSpline, bottomleft.x, bottomright.x,
                              intersectionPoints + SIDE_BOTTOM * 3, &intersectionPointsCount[SIDE_BOTTOM]);
    ac.splineLineIntersectionLat(A, B, leftSpline, topleft.y, bottomleft.y,
                              intersectionPoints + SIDE_LEFT * 3, &intersectionPointsCount[SIDE_LEFT]);
    ac.splineLineIntersectionLat(A, B, rightSpline, topright.y, bottomright.y,
                              intersectionPoints + SIDE_RIGHT * 3, &intersectionPointsCount[SIDE_RIGHT]);
    
    // now add all those points to array with sides!
    int intersectionPointsWithSidesCount = 0;
    SIntersectionPoint<double> intersectionPointsWithSides[14];
    for(int k = 0; k < 4; k++) {
        for(int l = 0; l < intersectionPointsCount[k]; l++) {
            if(out && count) {
                out[intersectionPointsWithSidesCount + *count].point = DHom(intersectionPoints[l + k * 3]);
                out[intersectionPointsWithSidesCount + *count].side = (CellSide)k;
                intersectionPointsWithSidesCount++;
            }
        }
    }
    
    // update counter
    if(count)*count += intersectionPointsWithSidesCount;
}


/// method to calculate area of given points with sides and connected splines
double VolumeCalculator::llPolygonWithSplinesAreaUnitSphere(SIntersectionPoint<double> *points, const int count, const double2 &center, const double *topSpline, const double *bottomSpline, const double *leftSpline, const double *rightSpline, const bool realWorldData) {
    
    // similiar to "linear" method, sort first points after angles relative to center point
#warning make here an assert check maybe for the center point! (must lie in polygon!)
    
#ifdef MATLAB
    double3 col = double3(0.05, 0.6, 0.05);
    double3 cols = double3(0.05, 0.05, 0.7);
    double3 colp = 0;
    ofs<<std::endl<<"scatter("<<center.x<<","<<center.y<<", 20, 'o', 'fill', 'MarkerFaceColor',["<<col.x<<" "<<col.y<<" "<<col.z<<"], 'MarkerEdgeColor', 'none'); hold on;\n";
#endif
    // make easy triangulation
    // no area for less than 3 points
    if(count < 3)return 0.0;
    
    // first step is to sort points
    // therefore go through points and calc their angle in relation to x axes! (because it can be easily calculated with atan2!)
    
    SAreaPoint *areaPoints = new SAreaPoint[count]; // holds angles between center-point and first point-center
    // calc angles using atan2 to get angles for 4 quadrants!
    for(int i = 0; i < count; i++) {
        double2 L = points[i].point - center;
        
        areaPoints[i].p = points[i];
        // atan2(y, x) !
        areaPoints[i].angle = atan2(L.y, L.x);
    }
    
    /*
    // sort with insertion sort
    for(int i = 1; i < count; i++) {
        for(int k = i; k > 0 && angles[k] < angles[k - 1]; k--) {
            swap(points[k], points[k - 1]);
            swap(angles[k], angles[k - 1]);
        }
    }
    */

    // sort after angles & sides!
    qsort(areaPoints, count, sizeof(SAreaPoint), compAreaPoints);
    
    int pointcount = count;

    double area = 0;
    double2 A = center;
    SIntersectionPoint<double> B;
    SIntersectionPoint<double> C;
    // now calc area based on triangulation
    // A remains fixed as
    for(int i = 0; i < pointcount; i++) {
        B = areaPoints[i].p;
        C = areaPoints[(i + 1) % pointcount].p;
        
#ifdef MATLAB
        mPoint(B.point, 30, colp);
        std::cout<<"point B: "<<B.point.toString()<<" angle:"<<rad2deg(areaPoints[i].angle)<<" @side: "<<B.side<<std::endl;
#endif
        
        // now the interesting part comes
        // iff B & C are on the same side use specialized function to integrate area of a triangle with a spline segment
        if(B.side == C.side && (B.side == SIDE_TOP || B.side == SIDE_BOTTOM || B.side == SIDE_LEFT || B.side == SIDE_RIGHT)) {
#ifdef MATLAB
            mLine(B.point, center, cols);
            mLine(C.point, center, cols);
#endif
            // use correct method for each side
            switch(B.side) {
                // funcs of lon
                case SIDE_TOP:
                    if(B.point.lon > C.point.lon)swap(B, C);
#ifdef MATLAB
                    ofs<<AreaCalculator::instance().printOutSplineCode(topSpline, B.point.lon, C.point.lon, 'x', cols);
                    ofs.flush();
#endif
                    area += llTriangleWithSplineLon(B.point, C.point, center, topSpline, realWorldData);
                    break;
                case SIDE_BOTTOM:
                    if(B.point.lon > C.point.lon)swap(B, C);
#ifdef MATLAB
                    ofs<<AreaCalculator::instance().printOutSplineCode(bottomSpline, B.point.lon, C.point.lon, 'x', cols);
                    ofs.flush();
#endif
                    area += llTriangleWithSplineLon(B.point, C.point, center, bottomSpline, realWorldData);
                    break;
                // funcs of lat
                case SIDE_LEFT:
                    if(B.point.lat > C.point.lat)swap(B, C);
#ifdef MATLAB
                    ofs<<AreaCalculator::instance().printOutSplineCode(leftSpline, B.point.lat, C.point.lat, 'y', cols, 2);
                    ofs.flush();
#endif
                    area += llTriangleWithSplineLat(B.point, C.point, center, leftSpline, realWorldData);
                    break;
                case SIDE_RIGHT:
                    if(B.point.lat > C.point.lat)swap(B, C);
#ifdef MATLAB
                    ofs<<AreaCalculator::instance().printOutSplineCode(rightSpline, B.point.lat, C.point.lat, 'y', cols);
                    ofs.flush();
#endif
                    area += llTriangleWithSplineLat(B.point, C.point, center, rightSpline, realWorldData);
                    break;
                default:
                    std::cout<<">> error: unknown side, something went really wrong!"<<std::endl;
                    break;
            }
        }
        else {
            // it is only a "simple" triangle --> use normal triangle function
            area += llTriangleUnitSphere(center, B.point, C.point, realWorldData);
            
#ifdef MATLAB
            mLine(B.point, C.point, col);
            mLine(center, B.point, col);
            mLine(center, C.point, col);
            ofs.flush();
#endif
        }
        
    }
    
    SafeDeleteA(areaPoints);
    
    assert(area >= 0);
    return area;
}

/// functions accepts values in degrees
double VolumeCalculator::splineLineIntegral(const double _a, const double _b, const double *spline, const double2& _A, const double2& _B, const bool realWorldData) {
    
    // convert to radians
    double a = deg2rad(_a);
    double b = deg2rad(_b);
    double2 A = vector2deg2rad(_A);
    double2 B = vector2deg2rad(_B);
    
    
    if(!realWorldData) {
        // not a necessary assertion...
        //assert(a < b);
        
        double aa = a * a;
        double bb = b * b;
        double aaa = a * aa;
        double bbb = b * bb;
        double aaaa = a * aaa;
        double bbbb = b * bbb;
        
        assert(!fequals(A.lon, B.lon));
        
        double m = (B.lat - A.lat) / (B.lon - A.lon);
        double c = A.lat - m * A.lon;
        
        // assert construction was ok
        assert(fequals(m * A.lon + c, A.lat));
        assert(fequals(m * B.lon + c, B.lat));
        
        // now integrate s(x)-g(x) (where s(x) is the spline's cubic polynomial!)
        double s[4] = {deg2rad(spline[0]), deg2rad(spline[1]), deg2rad(spline[2]), deg2rad(spline[3])};
        s[0] -= c;
        s[1] -= m;
        
        return (0.25 * s[3] * bbbb + 1/3.0 * s[2] * bbb + 0.5 * s[1] * bb + s[0] * b) -
        (0.25 * s[3] * aaaa + 1/3.0 * s[2] * aaa + 0.5 * s[1] * aa + s[0] * a);
    }
    else {
        assert(!fequals(A.lon, B.lon));
        
        double m = (B.lat - A.lat) / (B.lon - A.lon);
        double c = A.lat - m * A.lon;
        // use numerical integration

        // evaluate sin integral numerically via GSL!
        // set integrand
        double params[6] = {deg2rad(spline[0]), deg2rad(spline[1]), deg2rad(spline[2]), deg2rad(spline[3]), c, m};
        gsl_function F;
        F.function = &slf;
        F.params = &params;
        double inta = a;
        double intb = b;
        
        if(inta > intb)::swap(inta, intb);
        
        double error = 0;
        double area = 0;
        // use adaptive gauss integration method of GSL!
        gsl_integration_qags (&F, inta, intb, 0, EPSILON, 1000,
                              workspace, &area, &error);
        
        return area;
        
    }
}