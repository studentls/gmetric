//
//  MathFunctions.h
//  nemoExtractor
//
//  Created by Leonhard Spiegelberg on 09.08.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef nemoExtractor_MathFunctions_h
#define nemoExtractor_MathFunctions_h

#include "Vector2.h"
#include <complex>

/// interpolates in intveral x1 <-> x2 with a cubic spline, computing derivatives automatically
template<typename T, typename S> S cubicSpline(const T x, const T x0, const T x1, const T x2, const T x3, const S z0, const S z1, const S z2, const S z3){
    // check if x lies in interval
    assert(x >= x1);
    assert(x <= x2);
    
    // check if order is correct
    assert(x0 <= x1 && x1 <= x2 && x2 <= x3);
    
    //if(x == x1)return z1;
    //if(x == x2)return z2;
    
    // first step: compute derivatives at borders of interval
    T dx1 = 0.5 * ((z1 - z0) / (x1 - x0) + (z2 - z1) / (x2 - x1));
    T dx2 = 0.5 * ((z2 - z1) / (x2 - x1) + (z3 - z2) / (x3 - x2));
    
    // second step: calc Hermite functions after transforming x into [0, 1]
    T h = x2 - x1;
    
    T t = (x - x1) / h;
    T t2 = t * t;  // squared
    T t3 = t * t2; // cubic
    
    T H1 = 2 * t3 - 3*t2 + 1;
    T H2 = -2 * t3 + 3*t2;
    T H3 = t3 - 2*t2 + t;
    T H4 = t3 - t2;
    
    // third steo: calc actual value
    S res = z1 * H1 + z2 * H2 + h * dx1 * H3 + h * dx2 * H4;
    return res;
}

/// function constructs coefficients of a cubic splines' (local) polynomial
/// where p(t) = out[3] * t^3 + out[2] * t^2 + ... + out[0]
template<typename T> void cubicSplineCoefficients(const T x0, const T x1, const T x2, const T x3, const T z0, const T z1, const T z2, const T z3, T *out){
    
    // check for valid pointer
    assert(out);
    
    // check if order is correct
    assert(x0 <= x1 && x1 <= x2 && x2 <= x3);
    assert(x1 < x2);
    
    // first step: compute derivatives at borders of interval
    // if points are equal, set derivative to zero(only to get a value point!)
    T x2x1 = x2 - x1;
    T x3x2 = x3 - x2;
    T x1x0 = x1 - x0;
    T dx1 = x2x1 == 0 || x1x0 == 0 ? 0 : 0.5 * ((z1 - z0) / (x1x0) + (z2 - z1) / (x2x1));
    T dx2 = x3x2 == 0 || x2x1 == 0 ? 0 : 0.5 * ((z2 - z1) / (x2x1) + (z3 - z2) / (x3x2));
    
    // second step: calc Hermite functions after transforming x into [0, 1]
    T h = x2 - x1;
    
    assert(h != 0);
    
    // output coefficients
    T d = z1;
    T c = h * dx1;
    T b = -3.0 * z1 + 3 * z2 - 2 * h * dx1 - h * dx2;
    T a = 2 * z1 - 2 * z2 + h * dx1 + h * dx2;
    
    // apply transformation function inverse, and calc new coefficients!
    T aa = a / (h * h * h);
    T bb = b / (h * h);
    T cc = c / h;
    
    out[0] = -aa * x1 * x1 * x1 + bb * x1 * x1 - cc * x1 + d;
    out[1] = 3 * aa * x1 * x1 - 2 * bb * x1 + cc;
    out[2] = -3 * aa * x1 + bb;
    out[3] = aa;
    
}

/// horner schema to evaluate a polynomial
template<typename T> T evalPolynomial(const T t, const T *a, const unsigned int n) {
    // check pointer
    assert(a);
    
    T result = a[n - 1];
    for(int i=n-2; i >= 0 ; --i)
        result = result * t + a[i];
    return result;
}

/// implements Neville's Algorithm
template<typename T> T AitkenNeville(const T xq, T *x, T *y, const unsigned int n) {
    assert(n != 0);
    using namespace std;
    
    // copy initial values
    T *p = new T[n * n];
    
    // for debug reasons...
    memset(p, 0, sizeof(T) * n * n);
    
    for(int i = 0; i < n; i++)p[i] = y[i];
    
    for(int i = 0; i < n - 1; i++) {
        for(int j = 0; j < n - i - 1; j++) {
            // new version
            T p1 = p[j];
            T p2 = p[j + 1];
            T dp = p2 - p1;
            T xj = x[i + j + 1];
            T xi = x[j];
            T f = (xq - xi) / (xj- xi);
            T val = p1 + f * dp;
            
            p[j] = val;
        }
    }
    
    SafeDeleteA(p);
    
    return p[0];
}

/// implements Neville's Algorithm and destructs y array!
template<typename T> T AitkenNevilleInplace(const T xq, T *x, T *y, const unsigned int n) {
    assert(n != 0);
    
    for(int i = 0; i < n - 1; i++) {
        for(int j = 0; j < n - i - 1; j++) {
            // new version
            T p1 = y[j];
            T p2 = y[j + 1];
            T dp = p2 - p1;
            T xj = x[i + j + 1];
            T xi = x[j];
            T f = (xq - xi) / (xj- xi);
            T val = p1 + f * dp;
            
            y[j] = val;
        }
    }
    return y[0];
}


/// function to solve quadratic roots, for function x^2+px+q
template<typename T> std::complex<T> squareRootsPQ(T p, T q, std::complex<T> *out, int *numRoots) {
    
    assert(out);
    
    // calc discriminant
    T D = p*p - 4 * q;
    
    if(D < 0) {
        T p2 = p / 2.0;
        T sqrtexp = sqrt(q - p2 * p2);
        std::complex<T> x1 = std::complex<T>(-p2, +sqrtexp);
        std::complex<T> x2 = std::complex<T>(-p2, - sqrtexp);
        
        out[0] = x1;
        out[1] = x2;
        
        // two complex roots
        *numRoots = 2;
    }
    else {
        T p2 = p / 2.0;
        T sqrtexp = sqrt(p2 * p2 - q);
        std::complex<T> x1 = std::complex<T>(-p2 + sqrtexp, 0);
        std::complex<T> x2 = std::complex<T>(-p2 - sqrtexp, 0);
        
        out[0] = x1;
        out[1] = x2;
        
        if(D == 0) {
            // one(!) double real root
            *numRoots = 1;
        }
        else {
            // two real roots
            *numRoots = 2;
        }
    }
}

/// helper sign function
 template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


/// this function solves (real!) roots of a cubic polynomial after Cardano/Vieta
/// where the polynomial is p(t) = t^3+at^2+bt+c
inline int cubicRealRoots(double a, double b, double c, double *out, int *numRoots = NULL) {
    using namespace std;
    
    // assert output
    assert(out);
    int count = -1;
    
    double x[3]; // all roots
    
    // we assume that a, b, c, d are all real! no complex data!
    // also in general, a cubic polynomial has either 3 or 1 real root iff a,b,c,d are real!
    
    //algorithm after Numerical Recipes in C, pp. 228
    double Q = (a*a - 3.0  * b) / 9.0;
    double R = (2 * a* a * a - 9 * a * b + 27.0 * c) / 54.0;
    
    // three real roots?
    double Q3 = Q * Q * Q;
    if(R * R < Q3) {
        double theta = acos(R / sqrt(Q3));
        
        // roots
        double sqrtQ = sqrt(Q);
        x[0] = -2.0 * sqrtQ * cos(theta / 3.0) - a / 3.0;
        x[1] = -2.0 * sqrtQ * cos((theta + 2 * PI) / 3.0) - a / 3.0;
        x[2] = -2.0 * sqrtQ * cos((theta - 2 * PI) / 3.0) - a / 3.0;
        
        out[0] = x[0];
        out[1] = x[1];
        out[2] = x[2];
        count = 3;
        
    }
    // only one real root (note that complex roots are always pairwise in polynomials!)
    else {
        
        // compute A, B
        // note that R^2-Q^3 >= 0, due to if,else!
        double A = -sgn(R) * cbrt(abs(R) + sqrt(R * R - Q3));
        
        // calc B
        double B = A == 0 ? 0 : Q / A;
        
        x[0] = (A + B) - a / 3.0;
        
        out[0] = x[0];
        out[1] = NAN;
        out[2] = NAN;
        count = 1;
    }
    
    // out numroots if wished
    if(numRoots)*numRoots = count;
    return count;
}

/// this function solves (real!) roots of a cubic polynomial after Cardano/Vieta
/// where the polynomial is p(t) = at^3+bt^2+ct+d
inline int cubicRealRoots(double a, double b, double c, double d, double *out, int *numRoots = NULL) {
    
    int count = -1;
    // this function is for general a, b, c, d!
    // so check special cases!
    if(a == 0) {
        // quadratic function ?
        if(b == 0) {
            // linear function ?
            if(c == 0) {
                // const function
                out[0] = NAN;
                out[1] = NAN;
                out[2] = NAN;
                
                if(d == 0) {
                    out[0] = 0.0;
                    count = 1;
                }
                else count = 0;
            }
            else
            {
                // linear function
                out[0] = - d / c;
                out[1] = NAN;
                out[2] = NAN;
                count = 1;
            }
        }
        else {
            // quadratic function!
            // calc discriminant
            double D = c * c - 4.0 * b * d;
            if(D < 0){
                count = 0;
                out[0] = NAN;
                out[1] = NAN;
                out[2] = NAN;

            }else {
                // after numerical recipes in C, 3rd edition pp. 227
                double q = -0.5 * (c + sgn(c) * sqrt(D));
                out[0] = q / b;
                out[1] = d / q;
                out[2] = NAN;
                count = 2;
                if(D == 0) {
                    out[1] = NAN;
                    count--;
                }
            }
        }
    }
    else {
        // cubic function, solve!
        return cubicRealRoots(b / a, c / a, d / a, out, numRoots);
    }
    
    // output count if wished
    if(numRoots)*numRoots = count;
    return count;
}

/*
#error implement here general formula from wikipedia!!!
    using namespace std;
    
    assert(out);
    // special case a = 0; a,b = 0; a,b,c = 0 here!
    
    
    assert(a != 0);
    // first create depressed cubic polynomial
    // t^3+pt+q
    T p = (3 * a * c - b*b) / (3.0 * a * a);
    T q = (2*b*b*b - 9 * a * b * c + 27.0 * a*a*d) / (27.0 * a * a * a);
    
    // now apply vieta's substitution
    // t = w - p / (3w)
    // leading to equation
    // w^3 + q - p^3 / (27w^3) = 0
    // <=> w^6 + q * w^3 - p^3 / 27 = 0 (which is a sextic equation)
    // this quadratic formula can be solved
    
    // let be u = q, v = -p^3/27, solve (complex!) quadratic formula x^2+ux+v = 0 where x = w^3!
    T u = q;
    T v = - (p * p * p) / 27.0;
    
    complex<T> qrts[2]; // quadratic roots
    int numSquareRoots = 0;
    squareRootsPQ(u, v, qrts, &numSquareRoots);
    
    // now get w!
    T rcrt = cbrt(qrts[0].real()); // get cubic root of real part
    
    T phi = arg(qrts[0]);
    
    // get the three complex roots!
    complex<T> w1 = polar(rcrt, 1/3.0 * phi);
    complex<T> w2 = polar(rcrt, 1/3.0 * phi + 2.0/3.0 * PI);
    complex<T> w3 = polar(rcrt, 1/3.0 * phi - 2.0/3.0 * PI);
    
    // get roots of depressed polynomial
    complex<T> t1 = w1 - p / (3.0 * w1);
    complex<T> t2 = w2 - p / (3.0 * w2);
    complex<T> t3 = w3 - p / (3.0 * w3);
    
    // inverse Tschirnaus transformation yields roots of original polynomial
    complex<T> x[3];
    T fac = - b / (3 * a);
    x[0] = t1 + fac;
    x[1] = t2 + fac;
    x[2] = t3 + fac;
    
    *numRoots = 0;
    int pos = 0;
    // only output real roots
    for(int i = 0; i < 3; i++) {
        if(fequals(x[i].imag(), 0)) {
            out[pos++] = x[i].real();
        }
    }
    
    *numRoots = pos;
    // done!
}
*/

/// typedefs for double complex type
typedef std::complex<double> cdouble;
typedef std::complex<float> cfloat;
#endif
