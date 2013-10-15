//
//  TArray.h
//  gmetric
//
//  Created by Leonhard Spiegelberg on 20.08.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef TArray_h
#define TArray_h

#include "helper.h"


#define GM_CLAMPX 0x1
#define GM_CLAMPY 0x2
#define GM_MIRRORX 0x4
#define GM_MIRRORY 0x8
#define GM_CLAMPXY (GM_CLAMPX | GM_CLAMPY)
#define GM_CLAMP GM_CLAMPXY
#define GM_MIRRORXY (GM_MIRRORX | GM_MIRRORY)
#define GM_MIRROR GM_MIRRORXY

/// class to access gridded data in a 2D array over boundaries
template<typename T> class TArray {
public:
    //constructor
    TArray(const unsigned int _w, const unsigned int _h):w(_w), h(_h) {
        data = new T[w * h];
        mode = GM_CLAMP;
    }
    
    ~TArray() {
        if(data)SafeDeleteA(data);
    }
    
    // simple set & get based on NEMO indices!!!
    T get(const int _i, const int _j) {
        int i = _i;
        int j = _j;
        
        // change coords, due to mode
        if(mode & GM_CLAMPX)i = max(min(i, (int)w - 1), 0);
        if(mode & GM_CLAMPY)j = max(min(j, (int)h - 1), 0);
        if(mode & GM_MIRRORX) {
            while(i < 0)i += w;
            while(i > w - 1)i -= w;
        }
        if(mode & GM_MIRRORY) {
            while(j < 0)j += h;
            while(j > h - 1)j -= h;
        }
        
        assert(i + w * j < w * h && i + w * j >= 0);
        return data[i + w * j];
    }
    
    T get(const int _i) {
        int i = _i % w;
        int j = _i / w;
        return get(i, j);
    }
    
    void setMode(const int _mode) {mode = _mode;}
    
    void set(const int _i, const int _j, T val) {
        int i = _i;
        int j = _j;
        
        // change coords, due to mode
        if(mode & GM_CLAMPX)i = max(min(i, w - 1), 0);
        if(mode & GM_CLAMPY)j = max(min(j, h - 1), 0);
        if(mode & GM_MIRRORX) {
            while(i < 0)i += w;
            while(i >= w - 1)i -= w;
        }
        if(mode & GM_MIRRORY) {
            while(j < 0)j += h;
            while(j >= h - 1)j -= h;
        }
        
        assert(i + w * j < w * h && i + w * j >= 0);
        data[i + w * j] = val;
    }
    
    void set(const int i, T val) {
        assert(i < w * h && i >= 0);
        data[i] = val;
    }
    
    inline unsigned int width() {return w;}
    inline unsigned int height() {return h;}
    
    void transpose() {
        // transpose data
        arrayTranspose(data, w, h);
    }
    
    void swapdims() {
        swap(w, h);
    }
private:
    T *data;
    int w;
    int h;
    
    unsigned int mode;
};

#endif
