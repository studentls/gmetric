//
//  NCGrid.h
//  gmetric
//
//  Created by Leonhard Spiegelberg on 17.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#ifndef __gmetric__NCGrid__
#define __gmetric__NCGrid__

// netcdf
#include <netcdfcpp.h>
#include <netcdf.h>

// classes
#include "Vector2.h"
#include "Vector3.h"
#include "helper.h"
#include "MathFunctions.h"
#include "TArray.h"

// test case
//#include "NCGridTest.h"

// orientations of an array!
#define OM_ERROR        -1              // orientation error!
#define OM_NONE         0x0             // no orientation at all for coordinates!
#define OM_NEMO         0x1             // nemo orientation, y axis positive up, x axis positive right (elements get adressed via (x,y))
#define OM_NEMO90CW     0x2             // this is the standard orientation used by C++!
#define OM_NEMO180CW    0x3
#define OM_NEMO270CW    0x4
#define OM_NEMOMX       (0x8 | OM_NEMO) // nemo mirrored along x axis
#define OM_NEMO90CWMX   (0x8 | OM_NEMO90CW)
#define OM_NEMO180CWMX   (0x8 | OM_NEMO180CW)
#define OM_NEMO270CWMX   (0x8 | OM_NEMO270CW)

/// class to load a netcdf file with bathymetry, lon, lat and reorient data to default NEMO orientation
class NCGrid {
private:
    TArray<double2> *coordinates;
    TArray<double>  *bathymetry;
    int     orientation;
    
    /// helper function, to detect the orientation of data
    int     detectOrientation();
    
    /// functions to transform NEMO indices to indices of another array orientation and back vice versa
    int2 NEMOToOrientation(const int i, const int j, const int orientation);
    int2 OrientationToNEMO(const int i, const int j, const int orientation);
    
    /// rotate indices 90° CW(height is the original height)
    inline int2    rot90CW(int2 iv, const int height) {
        return int2(height - 1 - iv.y, iv.x);
    }
    
    /// rotate indices 90° CCW(width is the original width)
    inline int2    rot90CCW(int2 iv, const int width) {
        return int2(iv.y, width - 1 - iv.x);
    }
    
    /// rotate indices by 180° CW/CCW (width, height are original width/height)
    inline int2     rot180(int2 iv, const int width, const int height) {
        return int2(width - 1 - iv.x, height - 1 - iv.y);
    }
    
    /// mirror indices in X direction(width is the original width)
    inline int2    mirX(int2 iv, const int width) {
        return int2(width - 1 - iv.x, iv.y);
    }
    
    /// mirror indices in Y direction(height is the original height)
    inline int2    mirY(int2 iv, const int height) {
        return int2(iv.x, height -1 - iv.y);
    }
    
    // init arrays
    void    init(const int width, const int height);
    /// set data of index
    void setCoordinate(int i, int j, double2 val);
    void  setDepth(int i, int j, double val);
public:
    NCGrid() {
        coordinates = NULL;
        bathymetry = NULL;
        orientation = OM_NEMO; // per default NEMO orientation
    }
    
    ~NCGrid() {
        SafeDelete(coordinates);
        SafeDelete(bathymetry);
        orientation = OM_NONE;
    }
    
    /// open netcdf file and extract data by variable names
    int     open(std::string file,
                 std::string bathymetryVarName,
                 std::string lonVarName,
                 std::string        latVarName, const bool swapDimensions = false);
    
    /// get data of index
    double2 getCoordinate(int i, int j);
    double  getDepth(int i, int j);
    
    /// width(of array in NEMO orientation)
    int     width();
    
    /// height(of array in NEMO orientation)
    int     height();
    
    /// transform indices to original indices
    int2    getOriginalIndices(int i, int j);
    
    /// helper function to transpose internal data
    void    transpose();
    
    // test classes is nice to NCGrid, allow them private access :)
    friend class NCGridTest;
    friend class AreaCalculatorTest;
    friend class VolumeCalculatorTest;
};

#endif /* defined(__gmetric__NCGrid__) */
