//
//  Vector2.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 24.05.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef nemoExtractor_Vector2_h
#define nemoExtractor_Vector2_h

#include <cmath>
#include "helper.h"
#include "Vector3.h"
#include <vector>
#include <limits>

// forward declaration
template<typename T> class Vector3;

/// template class for 2D Vectors
template <typename T>
class Vector2 {
public:
    union {
        struct {
            T x,y;
        };
        struct {
            T lon, lat;
        };
    };
    
    Vector2():x(0), y(0)   {}
    Vector2(const T c):x(c), y(c) {}
    Vector2(const Vector2<T> & v):x(v.x), y(v.y){}
    Vector2(const Vector3<T> & v):x(v.x), y(v.y){}
    Vector2(const T _x, const T _y):x(_x), y(_y) {}
    
    
    inline T length() {
        return (T)sqrt((double)(x * x + y * y));
    }
    inline void			normalize() {
		T l = length();
		//secure in debug mode, avoid division by zero!
#ifdef _DEBUG
		if(l == 0.0)*this = Vector2<T>();
#endif
		x /= l; y /= l;
	}
    
    inline Vector2<T>		minus(const Vector2<T>& v)	{return Vector2<T>(x - v.x, y - v.y);}
    
	inline std::string toString() {
		std::stringstream s;
		s<<"("<<x<<","<<y<<")";
		return s.str();
	}
	inline bool equals(const Vector2<T>& v) {
		double eps = EPSILON;//0.0000001;
		Vector2<T> diff = this->minus(v);
		return (abs(diff.x) < eps && abs(diff.y) < eps);
	}
    
	inline bool operator == (const Vector2<T>& v) {
		double eps = EPSILON;//0.0000001;
		Vector2<T> diff = this->minus(v);
		return (abs(diff.x) < eps && abs(diff.y) < eps);
	}
    
    /// some more operators
	void operator =		(const Vector2<T>& v) {x =  v.x; y =  v.y; }
	void operator =		(const Vector3<T>& v) {x =  v.x; y =  v.y; }
	void operator +=	(const Vector2<T>& v)	{x += v.x; y += v.y; }
	void operator -=	(const Vector2<T>& v) {x -= v.x; y -= v.y; }
	void operator *=	(const double f)	{x *= f;   y *= f; }
	void operator /=	(const double f)	{x /= f;   y /= f; }
	void operator *=	(const float f)	{x *= f;   y *= f; }
	void operator /=	(const float f)	{x /= f;   y /= f; }
    void operator *=	(const int f)	{x *= f;   y *= f; }
	void operator /=	(const int f)	{x /= f;   y /= f; }
    
    // casting to Vector2
    operator Vector3<T>() {return Vector3<T>(x, y, 0);}
};

/// helper functions
template <typename T>
double vector2Length(const Vector2<T>& v) {
    return sqrt((double)(v.x * v.x + v.y * v.y));
}

template <typename T>
double vector2LengthSq(const Vector2<T>& v) {
    return v.x * v.x + v.y * v.y;
}

template <typename T>
T vector2Distance(const Vector2<T>& a, const Vector2<T>& b) {
    Vector2<T> v = a - b;
    return (T)sqrt((double)(v.x * v.x + v.y * v.y));
}

///operators
template <typename T> Vector2<T> operator + (const Vector2<T>& a, const Vector2<T>& b){return Vector2<T>(a.x + b.x, a.y + b.y);}
template <typename T> Vector2<T> operator - (const Vector2<T>& a, const Vector2<T>& b){return Vector2<T>(a.x - b.x, a.y - b.y);}
template <typename T> Vector2<T> operator - (const Vector2<T>& v)					  {return Vector2<T>(-v.x, -v.y);}
template <typename T> Vector2<T> operator * (const double d, const Vector2<T>& v)	  {return Vector2<T>(v.x * d, v.y * d);}
template <typename T> Vector2<T> operator * (const Vector2<T>& v, const double d)	  {return Vector2<T>(v.x * d, v.y * d);}
template <typename T> T operator * (const Vector2<T>& a, const Vector2<T>& b){return a.x * b.x + a.y * b.y;}
template <typename T> Vector2<T> operator / (const double d, const Vector2<T>& v)	  {return Vector2<T>(v.x / d, v.y / d);}
template <typename T> Vector2<T> operator / (const Vector2<T>& v, const double d)	  {return Vector2<T>(v.x / d, v.y / d);}
// for cout
template<typename T> std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
    Vector2<T> _v = v;
    os << _v.toString();
    return os;
}
    
/// typedef for GeoPoint
typedef Vector2<float> GeoPoint;

/// typedef to save writing time
typedef Vector2<float> float2;
typedef Vector2<double> double2;
typedef Vector2<int> int2;
typedef Vector2<unsigned int> uint2;

/////// functions

// min/max functions for vectors(component wise)
template<typename T> Vector2<T> vector2Min(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>(min(a.x, b.x), min(a.y, b.y));
}
template<typename T> Vector2<T> vector2Max(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>(max(a.x, b.x), max(a.y, b.y));
}

/// bilinear interpolator for a regular grid
template<typename T>	Vector2<T> vector2Lerp(const Vector2<T>& A, const Vector2<T>& B,
                                               const Vector2<T>& C, const Vector2<T>& D, const double x, const double y)
{
    // use two helper Points
    Vector2<T> P(A + x * (B - A));
    Vector2<T> Q(C + x * (D - C));
    
    // interpolate between helper Points
    return P + y * (Q - P);
}

/// nearest Neighbour Euclidean
/// A, B, C, D are the coordinates of the values a, b, c, d
/// x, y are the coordinates of the desired interpolation point
/// function will return value of point, which is closest to (x,y)
template<typename T, typename S>	S vector2NearestNeighbourEuclidean(const Vector2<T>& A, const Vector2<T>& B,
                                                                       const Vector2<T>& C, const Vector2<T>& D,
                                                                       const S a, const S b, const S c, const S d, const T x, const T y)
{
    // calc squared distances
    Vector2<T> P = Vector2<T>(x, y);
    double dAP = vector2LengthSq(A - P);
    double dBP = vector2LengthSq(B - P);
    double dCP = vector2LengthSq(C - P);
    double dDP = vector2LengthSq(D - P);
    
    // return corresponding values
    if(dAP > dBP) {
        if(dBP > dCP) {
            if(dCP > dDP)return d;
            else return c;
        }
        else {
            if(dBP > dDP)return d;
            else return b;
        }
    }
    else {
        if(dAP > dCP) {
            if(dCP > dDP)return d;
            else return c;
        }
        else {
            if(dAP > dDP)return d;
            else return a;
        }
    }
}

/// nearest Neighbour Manhattan distance based
/// A, B, C, D are the coordinates of the values a, b, c, d
/// x, y are the coordinates of the desired interpolation point
/// function will return value of point, which is closest to (x,y)
template<typename T, typename S>	S vector2NearestNeighbourManhattan(const Vector2<T>& A, const Vector2<T>& B,
                                                                       const Vector2<T>& C, const Vector2<T>& D,
                                                                       const S a, const S b, const S c, const S d, const T x, const T y)
{
    // calc squared distances
    Vector2<T> P = Vector2<T>(x, y);
    Vector2<T> AP = A - P;
    Vector2<T> BP = B - P;
    Vector2<T> CP = C - P;
    Vector2<T> DP = D - P;
    double dAP = abs(AP.x) + abs(AP.y);
    double dBP = abs(BP.x) + abs(BP.y);
    double dCP = abs(CP.x) + abs(CP.y);
    double dDP = abs(DP.x) + abs(DP.y);
    
    // return corresponding values
    if(dAP > dBP) {
        if(dBP > dCP) {
            if(dCP > dDP)return d;
            else return c;
        }
        else {
            if(dBP > dDP)return d;
            else return b;
        }
    }
    else {
        if(dAP > dCP) {
            if(dCP > dDP)return d;
            else return c;
        }
        else {
            if(dAP > dDP)return d;
            else return a;
        }
    }
}
/// bilinear Interpolation
/// x1, x2, y1, y2 are the coordinates of the values a, b, c, d in a regular grid
/// so a <-> (x1, y1), b <-> (x2, y1), c <-> (x1, y2), d <-> (x2, y2)
/// x, y are the coordinates of the desired interpolation point
/// function will return value of point, which is closest to (x,y)
template<typename T, typename S>	S bilinearInterpolation(const T x1, const T x2, const T y1, const T y2, const S a, const S b, const S c, const S d, const T x, const T y)
{
    // use formula from https://en.wikipedia.org/wiki/Bilinear_interpolation#Nonlinear
    
    // interpolation x-direction
    T ix2x1 = 1.0 / (x2 - x1);
    T p0 = (x2 - x) * ix2x1;
    T p1 = (x - x1) * ix2x1;
    
    S p = p0 * a + p1 * b;
    S q = p0 * c + p1 * d;
    
    // interpolate in y-direction
    T iy2y1 = 1.0 / (y2 - y1);
    return (y2 - y) * iy2y1 * p + (y - y1) * iy2y1 * q;
}


/// function to test if point is incident to line segment
template<typename T> bool pointOnLineSegment(const Vector2<T> A, const Vector2<T> B, const Vector2<T>& P) {
    Vector3<T> g = cross(Vector3<T>(A.x, A.y, 1), Vector3<T>(B.x, B.y, 1));
    
    // calculate dot product
    T d = dot(Vector3<T>(P.x, P.y, 1), g);
    
    if(abs(d) < EPSILON) {
        // check if P lies between A & B
        return partRatio(A, B, P) > -EPSILON; // part Ratio >= 0!
    }
    
    return false;
}

/// helper function to evaluate the line through (x1, y1) and (x2, y2) at position t!
template<typename T> T evalLine(const T x1, const T y1, const T x2, const T y2, const T x) {
    assert(!fequals(x1, x2));
    T a = (y2 - y1) / (x2 - x1);
    T b = y1  - a * x1;
    
    return a * x + b;
}

/// test for a cell, if a point P is contained or lies on boundary
template<typename T> bool pointInCell(const Vector2<T> *verts, const Vector2<T>& P) {
    
    // assert correct sorting!
    assert(verts);
    assert(verts[0].x < verts[1].x);    // correct NEMO orientation!
    assert(verts[3].x < verts[2].x);
    assert(verts[0].y > verts[3].y);
    assert(verts[1].y > verts[2].y);
    
    T leftVal = evalLine(verts[0].y, verts[0].x, verts[3].y, verts[3].x, P.lat);//evalPolynomial(p.lat, left, 4);
    // check if point lies within border (use <= if boundary shall be excluded from test!)
    if(P.lon < leftVal)return false;
    T rightVal = evalLine(verts[1].y, verts[1].x, verts[2].y, verts[2].x, P.lat);//evalPolynomial(p.lat, right, 4);
    if(P.lon > rightVal)return false;
    T topVal = evalLine(verts[0].x, verts[0].y, verts[1].x, verts[1].y, P.lon);//evalPolynomial(p.lon, top, 4);
    if(P.lat > topVal)return false;
    T bottomVal = evalLine(verts[3].x, verts[3].y, verts[2].x, verts[2].y, P.lon);
    if(P.lat < bottomVal)return false;
    
    return true;
}
    
/// returns barycentric coordinates for given triangle
template<typename T> Vector3<T> Barycentric(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const Vector2<T>& P) {
    Vector2<T> v0 = B - A;
    Vector2<T> v1 = C - A;
    Vector2<T> v2 = P - A;
    
    // dot product
    T d00 = v0 * v0;
    T d01 = v0 * v1;
    T d11 = v1 * v1;
    T d20 = v2 * v0;
    T d21 = v2 * v1;
    T invDenom = 1.0 / (d00 * d11 - d01 * d01);
    T v = (d11 * d20 - d01 * d21) * invDenom;
    T w = (d00 * d21 - d01 * d20) * invDenom;
    T u = 1 - v - w;
    
    return Vector3<T>(u, v, w);
}

/// test for a cell, if a point P is contained or lies on boundary of a triangle
template<typename T> bool pointInTriangle(const Vector2<T> *verts, const Vector2<T>& P) {
    // use barycentric coordinates for this task
    // triangle inside criteria is:
    // P inside triangle or on boundary <=> for all barycentric coords a_i: 0 <= a_i <= 1
    
    // original code from http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    // based on Cramer's rule
//    // Compute barycentric coordinates (u, v, w) for
//    // point p with respect to triangle (a, b, c)
//    void Barycentric(Point a, Point b, Point c, float &u, float &v, float &w)
//    {
//        Vector v0 = b - a, v1 = c - a, v2 = p - a;
//        float d00 = Dot(v0, v0);
//        float d01 = Dot(v0, v1);
//        float d11 = Dot(v1, v1);
//        float d20 = Dot(v2, v0);
//        float d21 = Dot(v2, v1);
//        float invDenom = 1.0 / (d00 * d11 - d01 * d01);
//        v = (d11 * d20 - d01 * d21) * invDenom;
//        w = (d00 * d21 - d01 * d20) * invDenom;
//        u = 1.0f - v - w;
//    }

//#error not yet implemented...
    Vector3<T> barycentric = Barycentric(verts[0], verts[1], verts[2], P);
    
    if(barycentric.x < 0)return false;
    if(barycentric.x > 1)return false;
    if(barycentric.y < 0)return false;
    if(barycentric.y > 1)return false;
    if(barycentric.z < 0)return false;
    if(barycentric.z > 1)return false;
    return true;
}

    
/// calc triangle area given by 3 points A, B, C
template<typename T> double triangleXYArea(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C) {
    
    // use that |cross(AB, AC)| = 2 times ABC!
    Vector3<T> AB = B - A;
    Vector3<T> AC = C - A;
    
    return 0.5 * vector3Length(cross(AB, AC));
}

/// haversine formula to calculate spherical distance in km between points given in lat/lon
template<typename T> double distanceLLdegOld(const Vector2<T>& A, const Vector2<T>& B, const double R = 1.0) {
    double halfdeltalon = deg2rad(B.lon - A.lon)  / 2.0;
    double halfdeltalat = deg2rad(B.lat - A.lat) / 2.0;
    
    double sqrtexpression = sin(halfdeltalat) * sin(halfdeltalat) + cos(deg2rad(A.lat)) * cos(deg2rad(B.lat)) * sin(halfdeltalon) * sin(halfdeltalon);
    
    // for numerical stability it is recommended to use atan2 but there is nor real difference, just ignore it!
    double angle = 2.0 * asin(sqrt(sqrtexpression));
    //double angle = 2.0 * atan2(sqrt(sqrtexpression), sqrt(1.0 - sqrtexpression));
    return angle * R;
}

/// updated version
template<typename T> double distanceLLdeg(const Vector2<T>& _A, const Vector2<T>& _B, const double R = 1.0) {
    Vector2<T> A = Vector2<T>(deg2rad(_A.lon), deg2rad(_A.lat));
    Vector2<T> B = Vector2<T>(deg2rad(_B.lon), deg2rad(_B.lat));
    
    
    double deltalon = B.lon - A.lon;
    
    double x = cos(B.lat) * sin(deltalon);
    double y = cos(A.lat) * sin(B.lat) - sin(A.lat) * cos(B.lat) * cos(deltalon);
    
    double denominator = sin(A.lat) * sin(B.lat) + cos(A.lat) * cos(B.lat) * cos(deltalon);
    
    double numerator = sqrt(x * x + y * y);
    
    double orthodrome = atan2(numerator, denominator);
    
    return R * orthodrome;
}

/// spherical triangle area in km^2
template<typename T> double triangleLLdegArea(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& C, const double R = 1.0) {
    // get side lengths on unit sphere
    double a = distanceLLdeg(B, C);
    double b = distanceLLdeg(A, C);
    double c = distanceLLdeg(A, B);
    
    double s = (a + b + c) / 2.0;
    
    // now calc spherical excess
    double sqrtexpression = sqrt(tan(s / 2.0) * tan((s - a) / 2.0) * tan((s - b) / 2.0) * tan((s - c) / 2.0));
    
    // excess
    double E = 4.0 * atan(sqrtexpression);
    
    // return area = E * R^2
    return E * R * R;
}

/// spherical triangle area in km^2
template<typename T> double triangleXYZArea(const Vector3<T>& A, const Vector3<T>& B, const Vector3<T>& C, const double R = 1.0) {
    // Maybe this could be also achieved a bit more faster, as angles could be angles between plane/great circles normals!
    Vector3<T> _A = vector3Normalize(A);
    Vector3<T> _B = vector3Normalize(B);
    Vector3<T> _C = vector3Normalize(C);
    
    
    // handle special case, if two points are equal
    if(_A == _B)return 0;
    if(_A == _C)return 0;
    if(_B == _C)return 0;
    
    // get side lengths on unit sphere
    double a = acos(_B * _C);
    double b = acos(_A * _C);
    double c = acos(_A * _B);
    
    double s = (a + b + c) / 2.0;
    
    // now calc spherical excess
    double sqrtexpression = sqrt(tan(s / 2.0) * tan((s - a) / 2.0) * tan((s - b) / 2.0) * tan((s - c) / 2.0));
    
    // excess
    double E = 4.0 * atan(sqrtexpression);
    
    // return area = E * R^2
    return E * R * R;
}

/// area of spherical convex polygon(sorts also points!)
template<typename T> double convexPolygonLLdegArea(Vector2<T> *points, const unsigned int count, const double R = 1.0) {
    // make easy triangulation
    
    if(count == 0)return 0.0;
    
    // only triangles & more supported! no "zweiecke"
    //assert(count >= 3);
    
    if(count == 3)
        return triangleLLdegArea(points[0], points[1], points[2], R);
    
    // first step is to sort points
    // therefore go through points and calc their angle in relation to y axes! (because it can be easily calculated with tan!)
    Vector2<T> center = Vector2<T>(0);
    for(int i = 0; i < count; i++) {
        center += points[i];
    }
    center = center / (double)count;
    
    // calculate angles
    double *angles = new double[count];
    for(int i = 0; i < count; i++) {
        angles[i] = atan2(points[i].lat - center.lat, points[i].lon - center.lon);
    }
    
    //    // print array
    //    std::cout<<" polygon points"<<std::endl<<"----"<<std::endl;
    //    for(int i = 0; i < count; i++) {
    //        std::cout<<"["<<i<<"] = "<<points[i].toString()<<" a="<<rad2deg(angles[i])<<std::endl;
    //    }
    
    // sort with insertion sort
    for(int i = 1; i < count; i++) {
        for(int k = i; k > 0 && angles[k] < angles[k - 1]; k--) {
            swap(points[k], points[k - 1]);
            swap(angles[k], angles[k - 1]);
        }
    }
    
    //    // print array
    //    std::cout<<" polygon points after sorting"<<std::endl<<"----"<<std::endl;
    //    for(int i = 0; i < count; i++) {
    //        std::cout<<"["<<i<<"] = "<<points[i].toString()<<" a="<<rad2deg(angles[i])<<std::endl;
    //    }
    
    int pointcount = count;
    // swap doublets to end
    /* for(int i = 1; i < count; i++){
     if(points[i] == points[i-1]) {
     // shift all elements by one index down
     for(int j = i; j < count - 1; j++) {
     points[j] = points[j +1];
     pointcount--;
     }
     }
     }*/
    
    //    // print array
    //    std::cout<<" polygon points after doublet removal"<<std::endl<<"----"<<std::endl;
    //    for(int i = 0; i < pointcount; i++) {
    //        std::cout<<"["<<i<<"] = "<<points[i].toString()<<" a="<<rad2deg(angles[i])<<std::endl;
    //    }
    
    
    
    double area = 0;
    Vector2<T> A = points[0];
    Vector2<T> B;
    Vector2<T> C;
    // now calc area based on triangulation
    for(int i = 1; i < pointcount - 1; i++) {
        B = points[i];
        C = points[i + 1];
        area += triangleLLdegArea(A, B, C, R);
    }
    
    SafeDeleteA(angles);
    
    return area;
}


/// area of spherical convex polygon(sorts also points!)
template<typename T> double convexPolygonXYZArea(Vector3<T> *points, const unsigned int count, const double R = 1.0) {
    // make easy triangulation
    
    if(count == 0)return 0.0;
    
    // only triangles & more supported! no "zweiecke"
    //assert(count >= 3);
    
    if(count == 3)
        return triangleXYZArea(points[0], points[1], points[2], R);
    
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
    
    //    // print array
    //    std::cout<<" polygon points"<<std::endl<<"----"<<std::endl;
    //    for(int i = 0; i < count; i++) {
    //        std::cout<<"["<<i<<"] = "<<points[i].toString()<<" a="<<rad2deg(angles[i])<<std::endl;
    //    }
    
    // sort with insertion sort
    for(int i = 1; i < count; i++) {
        for(int k = i; k > 0 && angles[k] < angles[k - 1]; k--) {
            swap(points[k], points[k - 1]);
            swap(angles[k], angles[k - 1]);
        }
    }
    
    //    // print array
    //    std::cout<<" polygon points after sorting"<<std::endl<<"----"<<std::endl;
    //    for(int i = 0; i < count; i++) {
    //        std::cout<<"["<<i<<"] = "<<points[i].toString()<<" a="<<rad2deg(angles[i])<<std::endl;
    //    }
    
    int pointcount = count;
    // swap doublets to end
    /* for(int i = 1; i < count; i++){
     if(points[i] == points[i-1]) {
     // shift all elements by one index down
     for(int j = i; j < count - 1; j++) {
     points[j] = points[j +1];
     pointcount--;
     }
     }
     }*/
    
    //    // print array
    //    std::cout<<" polygon points after doublet removal"<<std::endl<<"----"<<std::endl;
    //    for(int i = 0; i < pointcount; i++) {
    //        std::cout<<"["<<i<<"] = "<<points[i].toString()<<" a="<<rad2deg(angles[i])<<std::endl;
    //    }
    //
    
    
    double area = 0;
    Vector3<T> A = points[0];
    Vector3<T> B;
    Vector3<T> C;
    // now calc area based on triangulation
    for(int i = 1; i < pointcount - 1; i++) {
        B = points[i];
        C = points[i + 1];
        area += triangleXYZArea(A, B, C, R);
    }
    
    SafeDeleteA(angles);
    
    return area;
}

/// helper function to calculate part ratio FAST! (without check that points lie on same line!)
template<typename T> T partRatio(const Vector2<T>& A, const Vector2<T>& B, const Vector2<T>& P) {
    
    Vector2<T> AP = P - A;
    Vector2<T> PB = B - P;
    
    if(fequals(AP.x, 0) && fequals(AP.y, 0))return 0;
    if(fequals(PB.x, 0) && fequals(PB.y, 0))return INFINITY;
    
    // special case horizontal or vertical line
    if(!fequals(A.x, B.x))return AP.x / PB.x;
    else return AP.y / PB.y;
}

/// calculates intersection points between triangle and a quad
template<typename T> bool    quadrilateralTriangleIntersection(const Vector2<T> *abcPoints, const Vector2<T> *quadPoints, Vector2<T> *out, int *count) {
    
    bool intersection = false;
    int pos = 0;
    
    // homogenize
    Vector3<T> P[4];
    for(int i = 0; i < 4; i++)P[i] = Vector3<T>(quadPoints[i].x, quadPoints[i].y, 1);
    
    // calc homogenous points & triangle sides
    Vector3<T> s[3];
    Vector3<T> S[3];
    for(int i = 0; i < 3; i++)S[i] = Vector3<T>(abcPoints[i].x, abcPoints[i].y, 1);
    for(int i = 0; i < 3; i++)s[i] = cross(S[i], S[(i + 1) % 3]);
    
    
    // now calc intersection points
    for( int i = 0; i < 4; i++) {
        // line of quad side
        Vector3<T> h = cross(P[i], P[(i + 1) % 4]);
        
        for(int j = 0; j < 3; j++) {
            
            Vector3<T> HQ = cross(h, s[j]);
            
            // is Q inf point? => continue!
            if(fequals(HQ.z, 0))continue;
            
            Vector2<T> Q = Vector2<T>(HQ.x / HQ.z, HQ.y / HQ.z); // dehomogenize!
            
            // calc part ratios!
            // only fast method used!
            double prABC = partRatio(abcPoints[j], abcPoints[(j + 1) % 3], Q);
            double prQuad = partRatio(quadPoints[i], quadPoints[(i + 1) % 4], Q);
            
            // if both part ratios are positive, we have an intersection point on both triangle and quad segments!
            // note, that it is not allowed for them to be both zero! This would mean they fall together with an edge point!
            if(prABC >= 0 && prQuad >=  0 && (prABC > 0 || prQuad > 0)) {
                *(out + pos) = Vector2<T>(Q.x, Q.y);
                pos++;
                (*count)++;
                intersection = true;
            }
        }
    }
    
    return intersection;
}

/// calculates intersection points between triangle and a quad
template<typename T> bool    quadrilateralTriangleIntersection(const Vector2<T> *abcPoints, const Vector2<T> *quadPoints, Vector3<T> *out, int *count) {
    
    bool intersection = false;
    int pos = 0;
    
    // homogenize
    Vector3<T> P[4];
    for(int i = 0; i < 4; i++)P[i] = Vector3<T>(quadPoints[i].x, quadPoints[i].y, 1);
    
    // calc homogenous points & triangle sides
    Vector3<T> s[3];
    Vector3<T> S[3];
    for(int i = 0; i < 3; i++)S[i] = Vector3<T>(abcPoints[i].x, abcPoints[i].y, 1);
    for(int i = 0; i < 3; i++)s[i] = cross(S[i], S[(i + 1) % 3]);
    
    
    // now calc intersection points
    for( int i = 0; i < 4; i++) {
        // line of quad side
        Vector3<T> h = cross(P[i], P[(i + 1) % 4]);
        
        for(int j = 0; j < 3; j++) {
            
            Vector3<T> HQ = cross(h, s[j]);
            
            // is Q inf point? => continue!
            if(fequals(HQ.z, 0))continue;
            
            Vector2<T> Q = Vector2<T>(HQ.x / HQ.z, HQ.y / HQ.z); // dehomogenize!
            
            // calc part ratios!
            // only fast method used!
            double prABC = partRatio(abcPoints[j], abcPoints[(j + 1) % 3], Q);
            double prQuad = partRatio(quadPoints[i], quadPoints[(i + 1) % 4], Q);
            
            // if both part ratios are positive, we have an intersection point on both triangle and quad segments!
            // note, that it is not allowed for them to be both zero! This would mean they fall together with an edge point!
            if(prABC >= 0 && prQuad >=  0 && (prABC > 0 || prQuad > 0)) {
                *(out + pos) = HQ / HQ.z;
                pos++;
                (*count)++;
                intersection = true;
            }
        }
    }
    
    return intersection;
}

/// calculates intersection points between line and a quad,
/// line is given as element of RP^2
template<typename T> bool    quadrilateralLineIntersection(const Vector3<T> A, const Vector3<T> B, const Vector2<T> *quadPoints, Vector3<T> *out, int *count) {
    
    bool intersection = false;
    int pos = 0;
    
    // homogenize
    Vector3<T> P[4];
    for(int i = 0; i < 4; i++)P[i] = Vector3<T>(quadPoints[i].x, quadPoints[i].y, 1);
    
    // get line
    Vector3<T> line = cross(A, B);
    Vector2<T> _A = A / A.z;
    Vector2<T> _B = B / B.z;
    
    // now calc intersection points
    for( int i = 0; i < 4; i++) {
        // line of quad side
        Vector3<T> h = cross(P[i], P[(i + 1) % 4]);
        
        Vector3<T> HQ = cross(h, line);
        
        // is Q inf point? => continue!
        if(fequals(HQ.z, 0))continue;
        
        Vector2<T> Q = Vector2<T>(HQ.x / HQ.z, HQ.y / HQ.z); // dehomogenize!
        
        // calc part ratios!
        // only fast method used!
        double prAB = partRatio(_A, _B, Q);
        double prQuad = partRatio(quadPoints[i], quadPoints[(i + 1) % 4], Q);
        
        // if both part ratios are positive, we have an intersection point on both line and quad segments!
        // note, that it is not allowed for them to be both zero! This would mean they fall together with an edge point!
        if(prAB >= 0 && prQuad >=  0 && (prAB > 0 || prQuad > 0)) {
            *(out + pos) = HQ / HQ.z;
            pos++;
            (*count)++;
            intersection = true;
        }
    }
    
    return intersection;
}


/// calculates intersection points between triangle and a quad
template<typename T> bool    quadrilateralTriangleIntersectionXYZ(const Vector3<T> *abcPoints, const Vector3<T> *quadPoints, Vector3<T> *out, int *count) {
    
    bool intersection = false;
    int pos = 0;
    
    // calc homogenous triangle sides
    Vector3<T> s[3];
    Vector3<T> *P = quadPoints;
    Vector3<T> *S = abcPoints;
    for(int i = 0; i < 3; i++)s[i] = cross(S[i], S[(i + 1) % 3]);
    
    
    // now calc intersection points
    for( int i = 0; i < 4; i++) {
        // line of quad side
        Vector3<T> h = cross(P[i], P[(i + 1) % 4]);
        
        for(int j = 0; j < 3; j++) {
            
            Vector3<T> HQ = cross(h, s[j]);
            
            // is Q inf point? => continue!
            if(fequals(HQ.z, 0))continue;
            
            Vector2<T> Q = Vector2<T>(HQ.x / HQ.z, HQ.y / HQ.z); // dehomogenize!
            
            // calc part ratios!
            // only fast method used!
            double prABC = partRatio(abcPoints[j], abcPoints[(j + 1) % 3], Q);
            double prQuad = partRatio(quadPoints[i], quadPoints[(i + 1) % 4], Q);
            
            // if both part ratios are positive, we have an intersection point on both triangle and quad segments!
            // note, that it is not allowed for them to be both zero! This would mean they fall together with an edge point!
            if(prABC >= 0 && prQuad >=  0 && (prABC > 0 || prQuad > 0)) {
                *(out + pos) = HQ;
                pos++;
                (*count)++;
                intersection = true;
            }
        }
    }
    
    return intersection;
}


// converts vector deg2rad
template<typename T> Vector2<T> vector2deg2rad(const Vector2<T>& v) {
    return Vector2<T>(deg2rad(v.x), deg2rad(v.y));
}

template<typename T> Vector2<T> vector2rad2deg(const Vector2<T>& v) {
    return Vector2<T>(rad2deg(v.x), rad2deg(v.y));
}

/// calculates midpoint, given two lon/lat coords (source: http://gis.stackexchange.com/questions/55166/how-to-calculate-midpoint-given-north-south-latitude-and-east-west-longitude)
//Bx = cos(φ2).cos(Δλ)
//By = cos(φ2).sin(Δλ)
//φm = atan2( sin(φ1) + sin(φ2), √((cos(φ1)+Bx)² + By²) )
//λm = λ1 + atan2(By, cos(φ1)+Bx)
template<typename T> Vector2<T> midPointDeg(const Vector2<T>& _A, const Vector2<T>& _B) {
    Vector2<T> A = vector2deg2rad(_A);
    Vector2<T> B = vector2deg2rad(_B);
    
    T deltalon = B.lon - A.lon;
    T Bx = cos(B.lat) * cos(deltalon);
    T By = cos(B.lat) * sin(deltalon);
    
    T a = cos(A.lat) + Bx;
    T lat = atan2(sin(A.lat) + sin(B.lat), sqrt(a * a + By * By));
    T lon = A.lon + atan2(By, cos(A.lon) + Bx);
    
    return Vector2<T>(rad2deg(lon), rad2deg(lat));
}


/// little helper to assert spherical coords
template<typename T> void assertlonlat(const Vector2<T>& v) {
    assertlon(v.lon);
    assertlat(v.lat);
}

/// little helper to sort 3 points after longitude(so later A.lon < B.lon < C.lon!)
template<typename T> void sortABClon(Vector2<T>& A, Vector2<T>& B, Vector2<T>& C) {
    Vector2<T> h;
    if (A.lon > B.lon) {
        //swap(A, B);
        h = A;
        A = B;
        B = h;
    }
    if (B.lon > C.lon) {
        //Swap(B, C);
        h = B;
        B = C;
        C = h;
        
        if (A.lon > B.lon) {
            // Swap(A, B);
            h = A;
            A = B;
            B = h;
        }
    }
}


/// little helper to swap x & y coords
template<typename T> Vector2<T> swapXY(const Vector2<T> v) {
    return Vector2<T>(v.y, v.x);
}
    
/// helper to get rotated indices
template<typename T> Vector2<T> rotateIndices90CW(T i, T j, T w, T h) {
    return Vector2<T>(h - 1 - j, i);
}

/// helper to get rotated indices
template<typename T> Vector2<T> rotateIndices90CCW(T i, T j, T w, T h) {
    return Vector2<T>(j, w - 1 - i);
}

/// helper to get rotated indices
template<typename T> Vector2<T> mirrorIndicesX(T i, T j, T w, T h) {
    return Vector2<T>(w - 1 - i, j);
}
#endif
