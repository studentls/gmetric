#ifndef HELPER_HEADER_
#define HELPER_HEADER_

/// define some eps for floating point values
#define EPSILON 0.0000001

/// earth radius in km
#define EARTH_RADIUS 6378.137
#include <cassert>
#include <iostream>
#include <fstream>

// delete Helper
template <typename T>
void SafeDeleteA(T *array) {
    if(array)delete [] array;
    array = NULL;
}

template <typename T>
void SafeDelete(T *t) {
    if(t)delete t;
    t = NULL;
}

// helper struct for a bounding box
class BoundingBox {
public:
    double x;
    double y;
    double w;
    double h;
    
    BoundingBox():x(0), y(0), w(0), h(0) {}
    BoundingBox(const double _x,
                const double _y,
                const double _w,
                const double _h):x(_x),y(_y),w(_w),h(_h) {}
};

/// helper function to mirror array in x direction
template<typename T> void arrayMirrorX(T *a, const int w, const int h) {
    // in place!
    for(int i = 0; i < w / 2; i++) {
        for(int j = 0; j < h; j++) {
            // switch array elements
            int index1 = i + j * w;
            int index2 = (w - 1 - i) + j * w;
            
            // assert correct indices
            assert(index1 >= 0 && index1 < w * h);
            assert(index2 >= 0 && index2 < w * h);
            
            std::swap(a[index1], a[index2]);
        }
    }
}

/// helper function to mirror array in y direction
template<typename T> void arrayMirrorY(T *a, const int w, const int h) {
    // in place!
    for(int i = 0; i < w; i++) {
        for(int j = 0; j < h / 2; j++) {
            // switch array elements
            int index1 = i + j * w;
            int index2 = i + (h - 1 - j) * w;
            
            // assert correct indices
            assert(index1 >= 0 && index1 < w * h);
            assert(index2 >= 0 && index2 < w * h);
            
            std::swap(a[index1], a[index2]);
        }
    }
}

/// little helper to get gcd
template<typename T> T gcd(const T _a, const T _b)
{
    T c;
    T a = _a;
    T b = _b;
    while (a != 0) {
        c = a; a = b%a;  b = c;
    }
    return b;
}

/// m = w, n = h
/// helper function to transpose an array NOT IN PLACE!
/// note that after transpose w & h are swapped!
template<typename T> void arrayTranspose(T *a, int& w, int& h) {
    
    
    T *copy = new T[w * h];
    if(!copy)std::cout<<">> error: memory error during transpose!";
    else {
        // first copy
        for(int i = 0; i < w*h; i++)copy[i] = a[i];
        
        // transpose
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++) {
                a[i + j * w] = copy[j + i * h];
            }
        
        //delete copy
        SafeDeleteA(copy);
        
        std::swap(w, h);
    }
}

template<typename T> void arrayReverseElementsInCols(T *a, const int w, const int h) {
    
    // reverse all elements in every column...
    for(int i = 0; i < w; i++)
        for(int j = 0; j < h / 2; j++)std::swap(a[i + j * w], a[i + (h - 1 - j) * w]);
}

template<typename T> void arrayReverseElementsInRows(T *a, const int w, const int h) {
    
    // reverse all elements in every row...
    for(int j = 0; j < h; j++)
        for(int i = 0; i < w / 2; i++)std::swap(a[i + j * w], a[w - 1 - i + j * w]);
}

/// helper function to rotate an array clockwise
template<typename T> void arrayRotate90CW(T *a, int& w, int& h) {
    
    // (1) transpose
    arrayTranspose(a, w, h);
    
#warning for big arrays better use a more sophisticated version
    
    // (2) reverse elements in each row
    arrayReverseElementsInRows(a, w, h);
}

/// helper function to rotate an array counter clockwise
template<typename T> void arrayRotate90CCW(T *a, int& w, int& h) {
    
    // (1) transpose
    arrayTranspose(a, w, h);
    
    // (2) reverse elements in each column
    arrayReverseElementsInCols(a, w, h);
}

// helper to flip array
template <typename T>
void flipArray(T **a, int _w, int _h) {
    int wnew = _w;
    int hnew = _h;
    
    int w = _w;
    int h = _h;
    
    assert(*a);
    T *out = new T[hnew * wnew];
    assert(out);
    T *t = *a;
    for(int i = 0; i < w * h; i++) {
        int xold = i % w;
        int yold = i / w;
        int xnew = w - 1 - i % w;
        int ynew = h - 1 - i / w;
        
        out[xnew + ynew * wnew] = t[xold + yold * w];
        
    }
    delete [] (*a);
    *a = out;
}

// swap helper
template <typename T>
void swap(T& a, T& b) {
    T h = b;
    b = a;
    a = h;
}

// abs functions
template<typename T> T abs(const T& t) {
    return t > 0 ? t : -t;
}

// min/max functions
template <typename T>
T min(T a, T b) {
    return a < b ? a : b;
}

template <typename T>
T max(T a, T b) {
    return a > b ? a : b;
}

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170

#define PIDIV180 0.0174532925199432957692369076848861271344287188854172545609719144017100911460344944368224156963450948

// deg to rad converter
inline double deg2rad(const double deg) {
    //return deg / 180.0 * PI;
    return deg * PIDIV180;
}

// rad to deg
inline double rad2deg(const double rad) {
    return rad * 180.0 / PI;
}


// quickselect, used for median filter
// code from http://login2win.blogspot.ca/2011/06/quick-select.html

template<typename T> int partition(T* input, int p, int r)
{
    int pivot = input[r];
    
    while ( p < r )
    {
        while ( input[p] < pivot )
            p++;
        
        while ( input[r] > pivot )
            r--;
        
        if ( input[p] == input[r] )
            p++;
        else if ( p < r ) {
            int tmp = input[p];
            input[p] = input[r];
            input[r] = tmp;
        }
    }
    
    return r;
}

// p, r are left right and k is the kth rank
// r is the max element e.g.
template<typename T> T quick_select(T* input, int p, int r, int k)
{
    if ( p == r ) return input[p];
    int j = partition(input, p, r);
    int length = j - p + 1;
    if ( length == k ) return input[j];
    else if ( k < length ) return quick_select(input, p, j - 1, k);
    else  return quick_select(input, j + 1, r, k - length);
}


// median filter
template <typename T>
void medianFilter(T *a, const int w, const int h, const int kernel_size) {
    using namespace std;
    
    assert(a);
    
    // alloc mem for kernel
    T *kernel = new T[kernel_size * kernel_size];
    T *out = new T[w * h];
    
    // go through array
    for(int x = 0; x < w; x++)
    {
        for(int y = 0; y < h; y++) {
            
            // create kernel
            int pos = 0;
            for(int i = x - kernel_size / 2; i < x - kernel_size / 2 + kernel_size; i++)
                for(int j = y - kernel_size / 2; j < y - kernel_size / 2 + kernel_size; j++) {
                    int _i = ::min(w -1, ::max(0, i));
                    int _j = ::min(h - 1, ::max(0, h));
                    
                    assert(_i + _j * w < w * h);
                    kernel[pos++] = a[_i + _j * w];
                }
            // find median with quick select
            T val = quick_select(kernel, 0, kernel_size * kernel_size, kernel_size * kernel_size / 2);
            out[x + y * w] = val;
            
        }
        cout<<"   "<<(float)x * 100 /(float)w<<" %"<<endl;
    }
    delete [] a;
    a = out;
    
    delete [] kernel;
}

// a little helper to convert 16 bit from big endian to littl endian
inline short bigendian2littleendian(const short val) {
    short t = (val << 8) & 0xFF00;
    return t | ((val >> 8) & 0x00FF);
}

// helper to get only filename
inline std::string extractFilename(const std::string path) {
    std::string str = path;
    size_t index = str.rfind("/");
    str = str.substr(index + 1);
    return str;
}

inline bool fequals(float a, float b) {
    return ::abs(a - b) < EPSILON;
}

/// tridiagonal Matrix solver, destructive
/// the Matrix tries to solve
/// Ay = x
/// x contains first input and is used to store result
/// a is the subdiagonal and indexed from [1...N-1](so first element is dummy
/// b is the diagonal, indexed from [0....N-1]
/// c is the upperdiagonal, indexed from [1...N-2]
template <typename  T>
void solveTriDiagonalMatrix(T *x, const T *a, const T *b, const T *c, const int N) {
    // assertions
    assert(x);
    assert(a);
    assert(b);
    assert(c);
    assert(N >= 1);
    
    // we need dummy vector c'
    T *tc = new T[N];
    memcpy(tc, c, sizeof(T) * N);
    
    
    // start calculation
    tc[0] = c[0] / b[0];
    x[0]  = x[0] / b[0];
    
    // forward sweep
    for(int i = 1; i < N; i++) {
        T m = 1.0 / (b[i] - a[i] * tc[i - 1]);
        tc[i] = tc[i] * m;
        x[i] = (x[i] - a[i] * x[i - 1]) * m;
    }
    
    // backward substitution
    for(int i = N - 1; i-- > 0; ) {
        x[i] = x[i] - tc[i] * x[i + 1];
    }
    
    
    delete [] tc;
}

/// test if a value is contained in an array!
template<typename T> bool valueInArray(T *array, unsigned int count, T value) {
    for(int i = 0; i < count; i++) {
        if(array[i] == value) {
            return true;
        }
    }
    return false;
}

/// small helper to print array
template<typename T> void printArray(T *array, unsigned int count) {
    std::cout<<"|";
    for(int i = 0; i < count; i++)std::cout<<array[i]<<"|";
}

/// small helper to print squared array(numel = n * n) !
template<typename T> void printSquareArray(T *array, unsigned int n) {
    using namespace std;
    for(int j = 0; j < n; j++) {
        for(int i = 0; i < n; i++)cout<<array[i + j * n]<<"\t";
        cout<<endl;
    }
}

/// small helper to print rectangular array with parts
template<typename T> void printArray(T * array, unsigned int w, unsigned int h, int xstart = 0, int ystart = 0, int xend = - 1, int yend = -1) {
    using namespace std;
    if(xend < 0)xend = w - 1;
    if(yend < 0)yend = h - 1;
    
    assert(xstart >= 0);
    assert(ystart >= 0);
    assert(xend <= w - 1);
    assert(yend <= h - 1);
    
    for(int j = ystart; j <= yend; j++)
    {
        for(int i = xstart; i <= xend; i++)cout<<array[i + j * w]<<"\t";
        cout<<endl;
    }
}

/// small helper to print rectangular array with parts in NEMO style!
template<typename T> void printArrayNEMOStyle(T * array, unsigned int w, unsigned int h, int xstart = 0, int ystart = 0, int xend = - 1, int yend = -1) {
    using namespace std;
    if(xend < 0)xend = w - 1;
    if(yend < 0)yend = h - 1;
    
    assert(xstart >= 0);
    assert(ystart >= 0);
    assert(xend <= w - 1);
    assert(yend <= h - 1);
    
    assert(yend - ystart >= 1);
    assert(xend - xstart >= 1);
    for(int j = ystart; j <= yend; j++)
    {
        for(int i = xstart; i <= xend; i++) {
            // coord axes!
            if(j == yend && i == xstart)cout<<"^";
           // else if(j == yend && i == xstart)cout<<"|";
            else if(i == xstart || j == yend)cout<<" ";
            cout<<array[i + (h - 1 - j) * w]<<"\t";
        }
        cout<<endl;
    }
    cout<<"+ >"<<endl;
}

/// remove doublets in an array, maintain order!
/// pointer will be overwritten!
/// returns new count
template<typename T> unsigned int arrangeDoubletsToEnd(T *array, unsigned int count) {
    unsigned int doublets = 0;
    for(int j = 1; j < count - doublets; j++) {
        for(int i = 0; i < j; i++) {
            if(array[i] == array[j]) {
                // array at position j is the same as i! swap points!
                
                // find first place where a non doublet element is, at count-1-doublets element must not be doublet itself!
                bool found = false;
                while(!found && doublets != count - 2) {
                    // check if is doublet
                    bool isdoublet = false;
                    for(int k = 0; k <= j; k++) {
                        if(array[k] == array[count - 1 - doublets])isdoublet = true;
                    }
                    if(isdoublet)doublets++;
                    else found = true;
                }
                
                if(found) {
                    // swap elements
                    swap(array[j], array[count - 1 - doublets]);
                    doublets++;
                }
            }
        }
        //std::cout<<std::endl;
        //printArray(array, count);
    }
    return count - doublets;
}


/// little helper to assert lon/lat values lie in correct intervals
template<typename T> void assertlon(const T t) {
    assert(t  >= -PI);
    assert(t <= PI);
}
template<typename T> void assertlat(const T t) {
    assert(t  >= -PI/2.0);
    assert(t <= PI/2.0);
}

/// check if file exists
inline bool existsFile(const char* filename) {
    using namespace std;
    ifstream ifs(filename);
    return ifs;
}
#endif
