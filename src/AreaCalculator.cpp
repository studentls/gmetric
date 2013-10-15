//
//  AreaCalculator.cpp
//  gmetric
//
//  Created by Leonhard Spiegelberg on 18.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#include "AreaCalculator.h"

#define BUFFER_SIZE 14

/// function to implement the area metric
double AreaCalculator::linear(NCGrid *grid, const double2 _A, const double2 _B, const std::string filename, const bool realWorldData, const double R) {
    using namespace std;
    
//#error implement segmentfile
    // (1) first step is to get bounding box triangle in lon/lat coords
    double2 minABll = _A;
    double2 maxABll = _B;
    minABll = vector2Min(minABll, _B);
    maxABll = vector2Max(maxABll, _B);
    
    // (2) find bounding box indices!
    int si=-1,ei=-1,sj=-1,ej=-1;
    findBoundingBoxIndices(grid, minABll, maxABll, si, ei, sj, ej);
    
    // (3) now construct the points boundingBox
    //
    double2 pointsMin = grid->getCoordinate(si, sj);//_points[si + sj * w];
    double2 pointsMax = grid->getCoordinate(ei, ej);//_points[ei + ej * w];
    
    // assert this is a real bounding Box!
    assert(pointsMin.lon <= pointsMax.lon);
    assert(pointsMin.lat <= pointsMax.lat);
    
    // (4) now convert all lon/lat coordinates into (homogenized) lon/lat coords
    // first triangle coords
    double3 A = double3(_A.x, _A.y, 1);
    double3 B = double3(_B.x, _B.y, 1);
    
    // calc line between A & B
    double3 line = cross(A, B);
    
    // the area...
    double area = 0;
    
    // dehomogenized points in XYZ
    double2 lPoints[2];
    lPoints[0] = A / A.z;
    lPoints[1] = B / B.z;
    
    double2 p[4]; // edge points in XYZ
    
    //double3 *linePoints = new double3[BUFFER_SIZE];
    double3 linePoints[BUFFER_SIZE];
    int linePointsCount = 0;
    
    // if filename ist not empty, open output file and write comment line
    ofstream ofs;
    if(filename.length() != 0) {
        // does file exist alread? if so append!
        if(existsFile(filename.c_str())) {
            ofs.open(filename.c_str(), ios::app);
        }
        // print out header!
        else {
            ofs.open(filename.c_str());
            ofs<<"# i j lon1 lat1 lon2 lat2 distance area depth"<<endl;
        }
    }
    
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
            linePointsCount = 0;
            
            // now find out which points of A, B lie in Cell
            for(int k = 0; k < 2; k++) {
                if(pointInCell(p, lPoints[k]))linePoints[linePointsCount++] = double3(lPoints[k].x, lPoints[k].y, 1);
            }
            
            // now the most interesting part! determine intersection points and add to array
            quadrilateralLineIntersection(A, B, p, linePoints + linePointsCount, &linePointsCount);
            
            // if more than 2 points, sort!
            if(linePointsCount > 2)qsort((void*)linePoints, linePointsCount, sizeof(double3), compdouble3Points);
            
            // go through all line segments and add lengths!
            for(int k = 0; k < linePointsCount - 1; k++) {
                
                double2 lp1 = DHom(linePoints[k]);
                double2 lp2 = DHom(linePoints[k + 1]);
                
                // print warning if depth of zero is encountered!
                if(fequals(grid->getDepth(i, j), 0)) {
                    cout<<">> info: depth of zero encountered"<<endl;
                }
                
                // note that depth is in m^2! distance should return something in km!
                double distance = llDistanceUnitSphere(lp1, lp2, realWorldData);
                double segmentarea = distance * grid->getDepth(i, j);
                area += segmentarea;// * R * 0.001 later
                
                // write out segment if filename is valid
                if(filename.length() != 0) {
                    // get indices of original array!
                    int2 index = grid->getOriginalIndices(i, j);
                    
                    // just for conevenience, ensure that lp1.lon < lp2.lon!
                    if(lp1.lon > lp2.lon)::swap(lp1, lp2);
                    
                    // write out, format is
                    //# i j lon1 lat1 lon2 lat2 distance area depth
                    ofs<<index.x<<" "<<index.y<<" "<<lp1.lon<<" "<<lp1.lat<<" "<<lp2.lon<<" "<<lp2.lat<<" "<<(distance * R)<<" "<<(segmentarea * R * 0.001)<<" "<<grid->getDepth(i, j) * 0.001<<endl;
                }
                
            }
        }
    
    if(ofs.bad())cout<<" >> error: failure writing output segments file"<<endl;
    
    // return area scaled by R and 0.001 for the m to km conversion and unit sphere to earth sphere conversion!
    return area * R * 0.001;
}

/// function to return distance between two points given in lon/lat degree coordinates!
double AreaCalculator::llDistanceUnitSphere(const double2 A, const double2 B, bool realWorldData) {
    
    // convert to radians...
    double2 _A = vector2deg2rad(A);
    double2 _B = vector2deg2rad(B);
    
    // check if points are equal => do not solve integral then!
    if(_A == _B)return 0;
    
    if(!realWorldData) {
        // return ||.||_2 norm
        return (_A - _B).length();
    }
    else {
        // integrate
        double distance = 0;
        double error = 0;
        
        double3 line = cross(double3(_A.x, _A.y, 1), double3(_B.x, _B.y, 1));
        // set integrand
        double params[3] = {line.x, line.y, line.z};
        gsl_function F;
        F.function = &df;
        F.params = &params;
        double inta = abs(params[1]) < EPSILON ? _A.y : _A.x;
        double intb = abs(params[1]) < EPSILON ? _B.y : _B.x;
        
        if(inta > intb)::swap(inta, intb);
        
        // use adaptive gauss integration method of GSL!
        gsl_integration_qags (&F, inta, intb, 0, EPSILON, 1000,
                              workspace, &distance, &error);
        return distance;
    }
}

/// function to find bounding box indices
void AreaCalculator::findBoundingBoxIndices(NCGrid * grid, const double2 min, const double2 max, int& si, int & ei, int& sj, int & ej) {
    // (2) second step is to find indices of bounding box of points
    int w = grid->width();
    int h = grid->height();
    
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
        while(grid->getCoordinate(cursi + 1, j).lon < min.lon && cursi < w)cursi++;
        // find current ei
        while(grid->getCoordinate(curei - 1, j).lon > max.lon && curei >= 0)curei--;
        
        // global min/max along all columns!
        si = ::min(si, cursi);
        ei = ::max(ei, curei);
    }
    
    // now for lat
    for(int i = 0; i < w; i++)
    {
        int cursj = 0;
        int curej = h - 1;
        
        // find current si
        while(grid->getCoordinate(i, cursj + 1).lat < min.lat && cursj < w)cursj++;
        // find current ei
        while(grid->getCoordinate(i, curej - 1).lat > max.lat && curej >= 0)curej--;
        
        // global min/max along all columns!
        sj = ::min(sj, cursj);
        ej = ::max(ej, curej);
    }
}

/// function to implement the area metric using cubic splines as approximation for cells.
double AreaCalculator::spline(NCGrid *grid, const double2 _A, const double2 _B, const std::string filename, const bool realWorldData, const double R) {
    using namespace std;
    
    // (1) first step is to get bounding box triangle in lon/lat coords
    double2 minABll = _A;
    double2 maxABll = _B;
    minABll = vector2Min(minABll, _B);
    maxABll = vector2Max(maxABll, _B);
    
    // (2) find bounding box indices!
    int si=-1,ei=-1,sj=-1,ej=-1;
    findBoundingBoxIndices(grid, minABll, maxABll, si, ei, sj, ej);
    
    // (3) now construct the points boundingBox
    //
    double2 pointsMin = grid->getCoordinate(si, sj);//_points[si + sj * w];
    double2 pointsMax = grid->getCoordinate(ei, ej);//_points[ei + ej * w];
    
    // assert this is a real bounding Box!
    assert(pointsMin.lon <= pointsMax.lon);
    assert(pointsMin.lat <= pointsMax.lat);
    
    // (4) now convert all lon/lat coordinates into (homogenized) lon/lat coords
    // first triangle coords
    double3 A = double3(_A.x, _A.y, 1);
    double3 B = double3(_B.x, _B.y, 1);
    
    // calc line between A & B
    double3 line = cross(A, B);
    
    // the area...
    double area = 0;
    
    // dehomogenized points in XYZ
    double2 lPoints[2];
    lPoints[0] = A / A.z;
    lPoints[1] = B / B.z;
    
    double2 p[4]; // edge points in XYZ
    
    double3 *linePoints = new double3[BUFFER_SIZE];
    unsigned int linePointsCount = 0;
    
    // go through all cells (start with indices obtained in previous section)
    // as we access cells by i, i-1, j, j-1, shift indices by 1!
    si += 1; sj += 1; ei += 1; ej += 1;
    
    // ensure that there are enough points!
    assert(si < ei);
    assert(sj < ej);
    
    // if filename ist not empty, open output file and write comment line
    ofstream ofs;
    if(filename.length() != 0) {
        // does file exist alread? if so append!
        if(existsFile(filename.c_str())) {
            ofs.open(filename.c_str(), ios::app);
        }
        // print out header!
        else {
            ofs.open(filename.c_str());
            ofs<<"# i j lon1 lat1 lon2 lat2 distance area depth"<<endl;
        }
    }
    
#warning add here a better check for boundaries!
    
    // buffer for all intersection points
    SIntersectionPoint<double> intersectionPointsWithSides[14]; // maximum 14 needed!
    int intersectionPointsWithSidesCount = 0;
    
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
            double2 s0, s1, s2, s3;
            
            // !!! important, new selector for grid points should be implemented!
            
            
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
            
            
            // determine all line segments which intersect with cell sides
            vector<SLineSegment<double> > ls = splineSegmentsInsideCell(_A, _B, p[0], p[1], p[2], p[3],
                                                                        xTopSpline, xBottomSpline, yLeftSpline, yRightSpline);
            
            
            // calc distance
            double segmentarea = 0;
            double distance = 0;
            if(!ls.empty())
                for(vector<SLineSegment<double> >::const_iterator it = ls.begin(); it != ls.end(); ++it) {
                    
                    distance = llDistanceUnitSphere(it->X, it->Y, realWorldData);
                    
                    segmentarea = distance * grid->getDepth(i, j);
                    // note that depth is in m^2! distance should return something in km!
                    area += segmentarea;
                    
                    
                    // write out segment if filename is valid
                    if(filename.length() != 0) {
                        // get indices of original array!
                        int2 index = grid->getOriginalIndices(i, j);
                        
                        double2 lp1 = it->X;
                        double2 lp2 = it->Y;
                        
                        // just for conevenience, ensure that lp1.lon < lp2.lon!
                        if(lp1.lon > lp2.lon)::swap(lp1, lp2);
                        
                        // write out, format is
                        //# i j lon1 lat1 lon2 lat2 distance area depth
                        ofs<<index.x<<" "<<index.y<<" "<<lp1.lon<<" "<<lp1.lat<<" "<<lp2.lon<<" "<<lp2.lat<<" "<<(distance * R)<<" "<<(segmentarea * R * 0.001)<<" "<<grid->getDepth(i, j) * 0.001<<endl;
                    }
                    
                    // print warning if depth of zero is encountered!
                    if(fequals(grid->getDepth(i, j), 0)) {
                        cout<<">> info: depth of zero encountered"<<endl;
                    }
                }
        }
    
    SafeDeleteA(linePoints);
    
    // return area scaled by R and 0.001 for the m to km conversion and unit sphere to earth sphere conversion!
    return area * R * 0.001;
}

/// function to return all line segments, with value for splineArea intersection
std::vector<SLineSegment<double> > AreaCalculator::splineSegmentsInsideCell(const double2& A, const double2& B, const double2& topleft, const double2& topright, const double2& bottomright, const double2& bottomleft, double *xTopSpline, double *xBottomSpline, double *yLeftSpline, double *yRightSpline) {
    using namespace std;
    vector<SLineSegment<double> > res;
    
    // some asserts
    assert(xTopSpline);
    assert(xBottomSpline);
    assert(yLeftSpline);
    assert(yRightSpline);
    
    // intersect with splines
    // all in all,
    double3             intersectionPoints[4 * 3]; // 0-2 top, 3-5 bottom, 6-8 left, 9-11 right
    unsigned int        intersectionPointsCount[4] = {0, 0, 0, 0}; // 0 top, 1 bottom, 2 left, 3 right
    
    splineLineIntersectionLon(A, B, xTopSpline, topleft.x, topright.x,
                              intersectionPoints + SIDE_TOP * 3, &intersectionPointsCount[SIDE_TOP]);
    splineLineIntersectionLon(A, B, xBottomSpline, bottomleft.x, bottomright.x,
                              intersectionPoints + SIDE_BOTTOM * 3, &intersectionPointsCount[SIDE_BOTTOM]);
    splineLineIntersectionLat(A, B, yLeftSpline, topleft.y, bottomleft.y,
                              intersectionPoints + SIDE_LEFT * 3, &intersectionPointsCount[SIDE_LEFT]);
    splineLineIntersectionLat(A, B, yRightSpline, topright.y, bottomright.y,
                              intersectionPoints + SIDE_RIGHT * 3, &intersectionPointsCount[SIDE_RIGHT]);
    
    // now add all those points to array with sides!
    int intersectionPointsWithSidesCount = 0;
    SIntersectionPoint<double> intersectionPointsWithSides[14];
    for(int k = 0; k < 4; k++) {
        for(int l = 0; l < intersectionPointsCount[k]; l++) {
            intersectionPointsWithSides[intersectionPointsWithSidesCount].point = DHom(intersectionPoints[l + k * 3]);
            intersectionPointsWithSides[intersectionPointsWithSidesCount].side = (CellSide)k;
            intersectionPointsWithSidesCount++;
        }
    }
    
    // now find out which points of A, B lie in Cell
    if(llPointInSplineCell(A, yLeftSpline, yRightSpline, xTopSpline, xBottomSpline)) {
        intersectionPointsWithSides[intersectionPointsWithSidesCount].point = A;
        intersectionPointsWithSides[intersectionPointsWithSidesCount].side = SIDE_UNKNOWN;
        intersectionPointsWithSidesCount++;
    }
    if(llPointInSplineCell(B, yLeftSpline, yRightSpline, xTopSpline, xBottomSpline)) {
        intersectionPointsWithSides[intersectionPointsWithSidesCount].point = B;
        intersectionPointsWithSides[intersectionPointsWithSidesCount].side = SIDE_UNKNOWN;
        intersectionPointsWithSidesCount++;
    }
    
    // sort intersection points with sides array after lon(if lons are equal sort after lat!)
    qsort((void*)intersectionPointsWithSides, intersectionPointsWithSidesCount, sizeof(SIntersectionPoint<double>), compIntersectionPoints<double>);
    
    
    // calculate distance in cell!
    
    double distance = 0;
    SLineSegment<double> ls;
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
            double2 pointA = intersectionPointsWithSides[k].point;
            double2 pointB = intersectionPointsWithSides[k + 1].point;
            CellSide side = intersectionPointsWithSides[k].side;
            
            // now decide depending on sign of integral, if line between points lie in Cell or not!
            // x-based sides
            if(side == SIDE_TOP) {
                // top spline
                double intA = pointA.lon < pointB.lon ? pointA.lon : pointB.lon;
                double intB = pointA.lon < pointB.lon ? pointB.lon : pointA.lon;
                double integral = splineLineIntegral(intA, intB, xTopSpline, pointA, pointB);
                
                // if integral > 0 ==> add to distance!
                if(integral > 0) { //distance += Metrics::llDistance(pointA, pointB);
                    ls.X = pointA; ls.Y = pointB;
                    res.push_back(ls);
                }
                
            }
            if(side == SIDE_BOTTOM) {
                // bottom spline
                double intA = pointA.lon < pointB.lon ? pointA.lon : pointB.lon;
                double intB = pointA.lon < pointB.lon ? pointB.lon : pointA.lon;
                double integral = splineLineIntegral(intA, intB, xBottomSpline, pointA, pointB);
                
                // if integral < 0 ==> add to distance!
                if(integral < 0) { //distance += Metrics::llDistance(pointA, pointB);
                    ls.X = pointA; ls.Y = pointB;
                    res.push_back(ls);
                }
            }
            // y based sides
            if(side == SIDE_LEFT) {
                // left spline
                double intA = pointA.lat < pointB.lat ? pointA.lat : pointB.lat;
                double intB = pointA.lat < pointB.lat ? pointB.lat : pointA.lat;
                double integral = splineLineIntegral(intA, intB, yLeftSpline, swapXY(pointA), swapXY(pointB));
                
                // if integral > 0 ==> add to distance!
                if(integral > 0) { //distance += Metrics::llDistance(pointA, pointB);
                    ls.X = pointA; ls.Y = pointB;
                    res.push_back(ls);
                }
                
            }
            if(side == SIDE_RIGHT) {
                // right spline
                double intA = pointA.lat < pointB.lat ? pointA.lat : pointB.lat;
                double intB = pointA.lat < pointB.lat ? pointB.lat : pointA.lat;
                double integral = splineLineIntegral(intA, intB, yRightSpline, swapXY(pointA), swapXY(pointB));
                
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

/// helper function to solve integral between cubic polynomial of spline and line constructed from 2 points
/// a, b are the limits
double AreaCalculator::splineLineIntegral(const double a, const double b, const double *spline, const double2& A, const double2& B) {
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
        double s[4] = {spline[0], spline[1], spline[2], spline[3]};
        s[0] -= c;
        s[1] -= m;
        
        return (0.25 * s[3] * bbbb + 1/3.0 * s[2] * bbb + 0.5 * s[1] * bb + s[0] * b) -
        (0.25 * s[3] * aaaa + 1/3.0 * s[2] * aaa + 0.5 * s[1] * aa + s[0] * a);
}

/// helper which gives intersection point with a spline and a line
/// splineA, splineB define the interval for which the spline is defined
/// A, B define a line
bool    AreaCalculator::splineLineIntersectionLon(const double2 A, const double2 B, const double *spline, double splineALon, double splineBLon, double3 *out, unsigned int *count) {
    
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
            double3 IP = double3(A.lon, evalPolynomial(A.lon, spline, 4), 1);
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
        double a = (B.lat - A.lat) / (B.lon - A.lon);
        double b = A.lat - a * A.lon;
        
        // assert construction was ok
        assert(fequals(a * A.lon + b, A.lat));
        assert(fequals(a * B.lon + b, B.lat));
        
        // now find cubic roots of polynomial s(x)-g(x) (where s(x) is the spline's cubic polynomial!)
        double s[4] = {spline[0], spline[1], spline[2], spline[3]};
        s[0] -= b;
        s[1] -= a;
        
        // to prevent numerical errors, set values near to zero, to zero
        for(int i = 0; i < 4; i++)if(fequals(s[i], 0))s[i] = 0;
        
        double roots[3]; // there can be 3 roots!
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
                    double3 IP = double3(roots[i], a * roots[i] + b, 1);
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
bool    AreaCalculator::splineLineIntersectionLat(const double2& A, const double2& B, const double *spline, const double splineALat, const double splineBLat, double3 *out, unsigned int *count) {
    
    // save first, because lon/lats have to flipped back after computation
    double3      points[3]; // no more than 3 intersection points!
    unsigned int    pointscount = 0;
    bool res = splineLineIntersectionLon(double2(A.y, A.x), double2(B.y, B.x), spline, splineALat, splineBLat, points, &pointscount);
    
    // copy & flip coords, if no error occured
    if(res) {
        for(int i = 0; i < pointscount; i++) {
            // flip coordinates
            *(out + i) = double3(points[i].y, points[i].x, 1);
            (*count)++;
        }
    }
    
    return res;
}

std::string AreaCalculator::printOutSplineCode(const double *a, double start, double end, char mode, const double3& col, const double linewidth) {
    using namespace std;
    double h = 0.0025;
    stringstream xcords;
    stringstream ycords;
    
    //flip if necessary
    if(start > end)std::swap(start, end);
    
    for(int i = 0; i < (end - start) / h; i++) {
        // eval spline
        double t = start + i * h;
        double val = evalPolynomial(t, a, 4);
        
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
    double t = end;
    double val = evalPolynomial(t, a, 4);
    
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
    finalstr<<"plot(["<<xcords.str()<<"], ["<<ycords.str()<<"], 'Color', ["<<col.x<<" "<<col.y<<" "<<col.z<<"], 'LineWidth', "<<linewidth<<" ); hold on;\n";
    return finalstr.str();
}

/// helper to test, if point with lat, lon lies within cell bounded by splines top,left, right, bottom
/// point on boundary, returns true!!!
bool AreaCalculator::llPointInSplineCell(const double2& p, const double *left, const double * right, const double *top, const double *bottom) {
    // note: left & right are splines with lat as parameter
    //       top & bottom with lon as parameter
    
    double leftVal = evalPolynomial(p.lat, left, 4);
    // check if point lies within border (use <= if boundary shall be excluded from test!)
    if(p.lon < leftVal)return false;
    double rightVal = evalPolynomial(p.lat, right, 4);
    if(p.lon > rightVal)return false;
    double topVal = evalPolynomial(p.lon, top, 4);
    if(p.lat > topVal)return false;
    double bottomVal = evalPolynomial(p.lon, bottom, 4);
    if(p.lat < bottomVal)return false;
    
    return true;
}

bool AreaCalculator::pointsInAnyLinearCell(NCGrid *grid, const int si, const int sj, const int ei, const int ej, const double2& A, const double2& B) {
    bool foundA = false;
    bool foundB = false;
    double2 p[4];
    
    assert(si > 0);
    assert(sj > 0);
    assert(si <= ei);
    assert(sj <= ej);
    
    // loop through all points in BB
    for(int j=sj; j < ej; j++)
        for(int i=si; i < ei; i++)
        {
            // bottom layer can't be integrated
            if(i == 0 || j == 0)continue;
            
            // get Edge points
            p[0] = grid->getCoordinate(i - 1, j);
            p[1] = grid->getCoordinate(i,     j);
            p[2] = grid->getCoordinate(i,     j - 1);
            p[3] = grid->getCoordinate(i - 1, j - 1);
            
            if(pointInCell(p, A))foundA = true;
            if(pointInCell(p, B))foundB = true;
            
            // both found?
            if(foundA && foundB)return true;
        }
    return false;
}
