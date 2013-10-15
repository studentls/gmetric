//
//  Vector3.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 31.07.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef nemoExtractor_Vector3_h
#define nemoExtractor_Vector3_h

#include <cmath>
#include "helper.h"
#include "Vector2.h"
#include <sstream>

// forward declaration
template<typename T> class Vector2;

/// template class for 3D Vectors
template <typename T>
class Vector3 {
public:
    union {
        struct {
            T x,y,z;
        };
        struct {
            T lon, lat, r; // spherical coordinates
        };
    };
    
    Vector3():x(0), y(0), z(0)   {}
    Vector3(const T c):x(c), y(c), z(c) {}
    Vector3(const Vector3<T> & v):x(v.x), y(v.y), z(v.z){}
    Vector3(const Vector2<T> & v):x(v.x), y(v.y), z(0){}
    Vector3(const T _x, const T _y):x(_x), y(_y), z(0) {}
    Vector3(const T _x, const T _y, const T _z):x(_x), y(_y), z(_z) {}
    
    
    inline T length() {
        return (T)sqrt((double)(x * x + y * y + z * z));
    }
    inline void			normalize() {
		T l = length();
		//secure in debug mode, avoid division by zero!
#ifdef _DEBUG
		if(l == 0.0)*this = Vector3<T>();
#endif
		x /= l; y /= l;
	}
    
    inline Vector3<T>		minus(const Vector3<T>& v)	{return Vector3<T>(x - v.x, y - v.y, z - v.z);}
    
    inline Vector3<T>       cross(const Vector3<T>& v) {return Vector3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);}
    
	inline std::string toString() {
		std::stringstream s;
		s<<"("<<x<<","<<y<<","<<z<<")";
		return s.str();
	}
	inline bool equals(const Vector3<T>& v) {
		double eps = EPSILON;//0.0000001;
		Vector3<T> diff = this->minus(v);
		return (abs(diff.x) < eps && abs(diff.y) < eps);
	}
    
	inline bool operator == (const Vector3<T>& v) {
		double eps = EPSILON;//0.0000001;
		Vector3<T> diff = this->minus(v);
		return (abs(diff.x) < eps && abs(diff.y) < eps);
	}
    
    /// some more operators
	void operator =		(const Vector2<T>& v) {x =  v.x; y =  v.y; z = 0;}
	void operator =		(const Vector3<T>& v) {x =  v.x; y =  v.y; z = v.z;}
	void operator +=	(const Vector3<T>& v)	{x += v.x; y += v.y; z += v.z;}
	void operator -=	(const Vector3<T>& v) {x -= v.x; y -= v.y; z -= v.z;}
	void operator *=	(const double f)	{x *= f;   y *= f; z *= f;}
	void operator /=	(const double f)	{x /= f;   y /= f; z /= f;}
	void operator *=	(const float f)	{x *= f;   y *= f; z *= f;}
	void operator /=	(const float f)	{x /= f;   y /= f; z /= f;}
    void operator *=	(const int f)	{x *= f;   y *= f; z *= f;}
	void operator /=	(const int f)	{x /= f;   y /= f; z /= f;}
    
    // casting to Vector2
    operator Vector2<T>() {return Vector2<T>(x, y);}
    
};

/// helper functions
template <typename T>
double vector3Length(const Vector3<T>& v) {
    return sqrt((double)(v.x * v.x + v.y * v.y + v.z * v.z));
}

template <typename T>
double vector3LengthSq(const Vector3<T>& v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

template <typename T>
T vector3Distance(const Vector3<T>& a, const Vector3<T>& b) {
    Vector3<T> v = a - b;
    return (T)sqrt((double)(v.x * v.x + v.y * v.y + v.z * v.z));
}

template<typename T> Vector3<T> vector3Normalize(const Vector3<T>& v) {
    double length = vector3Length(v);
    
    assert(length != 0.0);
    return v / length;
}
// alternative, cross function directly
template<typename T> Vector3<T> cross(const Vector3<T>& a, const Vector3<T>& b) {
    return Vector3<T>(a.y * b.z - a.z * b.y,
                      a.z * b.x - a.x * b.z,
                      a.x * b.y - a.y * b.x);
}

// dot product
template<typename T> T dot(const Vector3<T>& a, const Vector3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// min/max functions for vectors(component wise)
template<typename T> Vector3<T> vector3Min(const Vector3<T>& a, const Vector3<T>& b) {
    return Vector3<T>(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}
template<typename T> Vector3<T> vector3Max(const Vector3<T>& a, const Vector3<T>& b) {
    return Vector3<T>(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}

///operators
template <typename T> Vector3<T> operator + (const Vector3<T>& a, const Vector3<T>& b){return Vector3<T>(a.x + b.x, a.y + b.y, a.z + b.z);}
template <typename T> Vector3<T> operator - (const Vector3<T>& a, const Vector3<T>& b){return Vector3<T>(a.x - b.x, a.y - b.y, a.z - b.z);}
template <typename T> Vector3<T> operator - (const Vector3<T>& v)					  {return Vector3<T>(-v.x, -v.y, -v.z);}
template <typename T> Vector3<T> operator * (const double d, const Vector3<T>& v)	  {return Vector3<T>(v.x * d, v.y * d, v.z * d);}
template <typename T> Vector3<T> operator * (const Vector3<T>& v, const double d)	  {return Vector3<T>(v.x * d, v.y * d, v.z * d);}
template <typename T> T operator * (const Vector3<T>& a, const Vector3<T>& b){return a.x * b.x + a.y * b.y + a.z * b.z;}
template <typename T> Vector3<T> operator / (const double d, const Vector3<T>& v)	  {return Vector3<T>(v.x / d, v.y / d, v.z / d);}
template <typename T> Vector3<T> operator / (const Vector3<T>& v, const double d)	  {return Vector3<T>(v.x / d, v.y / d, v.z / d);}
// for cout
template<typename T> std::ostream& operator<<(std::ostream& os, const Vector3<T>& v)
{
    Vector3<T> _v = v;
    os << _v.toString();
    return os;
}



/// typedef to save writing time
typedef Vector3<float> float3;
typedef Vector3<double> double3;



/// helper functions

/// convert spherical to xyz coords (Attention!!! uses not phyiscal standard, but instead more geographical interpretation!
template<typename T> Vector3<T> spherical2xyz(const Vector2<T> p) {
    assert(p.lat < PI/2 && p.lat > -PI/2 );
    assert(p.lon < PI && p.lon > -PI);
    
    return Vector3<T>(cos(p.lon) * cos(p.lat),sin(p.lon) * cos(p.lat),sin(p.lat));
}

/// convert xyz to spherical lon/lat (Attention!!! xyz vector has to be a unit vector!
template<typename T> Vector2<T> xyz2spherical(const Vector3<T> p) {
    return Vector2<T>(atan2(p.y, p.x),asin(p.z));
}

/// normalize vector
template<typename T> Vector3<T> vector3Normaliye(const Vector3<T> v) {
    T length = vector3Length(v);
    assert(fequals(0, length));
    return v / length;
    
}
///// convert spherical to xyz coords after physical standard
//template<typename T> Vector3<T> spherical2xyz(const Vector2<T> p) {
//    assert(p.lat < PI);
//    assert(p.lon < PI);
//    
//    return Vector3<T>(cos(p.lon) * sin(p.lat), sin(p.lon) * sin(p.lat), cos(p.lat));
//}
//
///// convert xyz to spherical lon/lat (Attention!!! xyz vector has to be a unit vector!
//template<typename T> Vector2<T> xyz2spherical(const Vector3<T> p) {
//    return Vector2<T>(atan2(p.y, p.x),acos(p.z));
//}

/// new function, inspired through https://groups.google.com/forum/#!topic/de.sci.mathematik/CLfwGm2hDwo
/// where spherical excess is
/// E = 2*atan [(a x b)*c)/(1+a*b+a*c+b*c)]
template<typename T> T xyzTriangleArea(const Vector3<T>& _A, const Vector3<T>& _B, const Vector3<T>& _C, const T R = 1.0) {
    // normalize
    Vector3<T> A = vector3Normalize(_A);
    Vector3<T> B = vector3Normalize(_B);
    Vector3<T> C = vector3Normalize(_C);
    
    T spat = cross(A, B) * C;
    T ab = A * B;
    T ac = A * C;
    T bc = B * C;
    T E = 2 * atan2(spat, 1 + ab + ac + bc);
    
    return E * R * R;
}

/// great circle distance
template<typename T> T xyzDistance(const Vector3<T>& _A, const Vector3<T>& _B, const T R = 1.0) {
    // normalize
    Vector3<T> A = vector3Normalize(_A);
    Vector3<T> B = vector3Normalize(_B);
    
    return (A * B) * R;
}

/// little helper to sort 3 points after longitude(so later A.lon < B.lon < C.lon!)
template<typename T> void sortABClon(Vector3<T>& A, Vector3<T>& B, Vector3<T>& C) {
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

/// little helper to homogenize coords
template<typename T> Vector3<T> Hom(const Vector2<T>& v) {
    return Vector3<T>(v.x, v.y, 1);
}
/// little helper to Dehomogenize coords
template<typename T> Vector2<T> DHom(const Vector3<T>& v) {
    //assert(v.z != 0);
    return v.z == 0 ? Vector2<T>(v.x, v.y)  : Vector2<T>(v.x / v.z, v.y / v.z);
}
#endif
