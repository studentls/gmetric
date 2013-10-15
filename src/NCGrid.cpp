//
//  NCGrid.cpp
//  gmetric
//
//  Created by Leonhard Spiegelberg on 17.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#include "NCGrid.h"

/// open netCDF file and pass data to array structure
/// returns positive value if no error occured
int NCGrid::open(std::string file, std::string bathymetryVarName, std::string lonVarName, std::string latVarName, const bool swapDimensions) {
    using namespace std;
    
    // load natl file, to extract positions
    NcFile ncFile = NcFile(file.c_str(), NcFile::ReadOnly);
    
    NcVar *lonVar = ncFile.get_var(lonVarName.c_str());
    NcVar *latVar = ncFile.get_var(latVarName.c_str());
    NcVar *depthVar = ncFile.get_var(bathymetryVarName.c_str());
    
    // first secure, that there are two dimensions
    if(lonVar->num_dims() != 2 || latVar->num_dims() != 2 || depthVar->num_dims() != 2) {
        cout<<">> error: dimension error, variables "<<bathymetryVarName<<","<<lonVarName<<","<<latVarName<<" have to be two dimensional arrays"<<endl;
        return -1;
    }
    
    // retrieve dimensions from vars
    NcDim *lonDims[2];
    NcDim *latDims[2];
    NcDim *depthDims[2];
    
    for(int i = 0; i < 2; i++) {
        lonDims[i] = lonVar->get_dim(i);
        latDims[i] = latVar->get_dim(i);
        depthDims[i] = depthVar->get_dim(i);
    }
    // print out dimensions:
    cout<<">> lon:        '"<<lonDims[0]->name()<<"' x '"<<lonDims[1]->name()<<"'\t"<<lonDims[0]->size()<<" x "<<lonDims[1]->size()<<endl;
    cout<<">> lat:        '"<<latDims[0]->name()<<"' x '"<<latDims[1]->name()<<"'\t"<<latDims[0]->size()<<" x "<<latDims[1]->size()<<endl;
    cout<<">> Bathymetry: '"<<depthDims[0]->name()<<"' x '"<<depthDims[1]->name()<<"'\t"<<depthDims[0]->size()<<" x "<<depthDims[1]->size()<<endl;
    
    // assert, that arrays have compatible dimensions!
    for(int i = 0; i < 2; i++)
        if(!(lonDims[i]->size() == latDims[i]->size() && lonDims[i]->size() == depthDims[i]->size())) {
            cout<<">> error: dimension error, variable dimension do not fit"<<endl;
            return -3;
        }
    
    // beyond this point, dimensions of variables do agree!
    int width = (int)lonDims[0]->size();
    int height = (int)lonDims[1]->size();
    
    // init internal arrays
    init(width, height);
    
    // regarding the type of the netcdf variable extract data
    if(lonVar->type() == ncFloat) {
        float *tmp = new float[width * height];
        lonVar->get(tmp, width, height);
        for(int i = 0; i < width * height; i++) {
            double2 tmpcoord = coordinates->get(i);
            tmpcoord.lon = (double)tmp[i];
            coordinates->set(i, tmpcoord);
        }
        SafeDeleteA(tmp);
    } else if(lonVar->type() == ncDouble) {
        double *tmp = new double[width * height];
        lonVar->get(tmp, width, height);
        for(int i = 0; i < width * height; i++) {
            double2 tmpcoord = coordinates->get(i);
            tmpcoord.lon = (double)tmp[i];
            coordinates->set(i, tmpcoord);
        }
        SafeDeleteA(tmp);
    }
    else {
        cout<<">> error: type error for variable "<<lonVar->name()<<" only supported type for the input grid are ncFloat or ncDouble"<<endl;
        return -2;
    }
    
    if(latVar->type() == ncFloat) {
        float *tmp = new float[width * height];
        latVar->get(tmp, width, height);
        for(int i = 0; i < width * height; i++) {
            double2 tmpcoord = coordinates->get(i);
            tmpcoord.lat = (double)tmp[i];
            coordinates->set(i, tmpcoord);
        }
        SafeDeleteA(tmp);
    } else if(latVar->type() == ncDouble) {
        double *tmp = new double[width * height];
        latVar->get(tmp, width, height);
        for(int i = 0; i < width * height; i++) {
            double2 tmpcoord = coordinates->get(i);
            tmpcoord.lat = (double)tmp[i];
            coordinates->set(i, tmpcoord);
        }
        SafeDeleteA(tmp);
    }
    else {
        cout<<">> error: type error for variable "<<latVar->name()<<" only supported type for the input grid are ncFloat or ncDouble"<<endl;
        return -2;
    }
    
    if(depthVar->type() == ncFloat) {
        float *tmp = new float[width * height];
        depthVar->get(tmp, width, height);
        for(int i = 0; i < width * height; i++)bathymetry->set(i, (double)tmp[i]);
        SafeDeleteA(tmp);
    } else if(depthVar->type() == ncDouble) {
        double *tmp = new double[width * height];
        depthVar->get(tmp, width, height);
        for(int i = 0; i < width * height; i++)bathymetry->set(i, (double)tmp[i]);
        SafeDeleteA(tmp);
    }
    else {
        cout<<">> error: type error for variable "<<depthVar->name()<<" only supported type for the input grid are ncFloat or ncDouble"<<endl;
        return -2;
    }
    
    // update all coords, such that coords are in range!
    for(int i = 0; i < width * height; i++) {
        double2 coord = coordinates->get(i);
        if(coord.x > 180)coord.x -= 360;
        if(coord.x < -180)coord.x += 360;
        if(coord.y > 90 || coord.y < -90) {
            cout<<">> error: invalid latitudes found"<<endl;
            return -19;
        }
        
        // write back
        coordinates->set(i, coord);
    }
    
    // swap coordinate dimensions if specified...
    if(swapDimensions) {
        coordinates->swapdims();
        bathymetry->swapdims();
    }
    
    cout<<">> info: loaded grid of "<<this->width()<<"x"<<this->height()<<" points"<<endl;
    
    // now comes the orientation code...
    orientation = detectOrientation();
    
    // error in coordinates?
    if(orientation == OM_NONE) {
        cout<<">> error: coordinate error - data has no orientation"<<endl;
        return -5;
    }
    
    // there should be no -1 returned indicating failure of routine!
    assert(orientation >= 0);
    
    return 1;
}

/// detects orientation of the input data and errors by scanning coordinates
int NCGrid::detectOrientation() {
    
    using namespace std;
    // check if x-axis is ascending/descending
    bool errX = false;
    bool xLonAsc = false;
    bool xLonDesc = false;
    bool xLatAsc = false;
    bool xLatDesc = false;
    for(int j = 0; j < coordinates->height(); j++) { // go through all y layers
        for(int i = 0; i < coordinates->width() - 1; i++) {
            double2 p1 = coordinates->get(i,     j);
            double2 p2 = coordinates->get(i + 1, j);
            if(p1.lon < p2.lon)xLonAsc = true;
            if(p1.lon > p2.lon)xLonDesc = true;
            if(p1.lat < p2.lat)xLatAsc = true;
            if(p1.lat > p2.lat)xLatDesc = true;
        }
    }
    
    // check if y-axis is ascending/descending
    bool errY = false;
    bool yLonAsc = false;
    bool yLonDesc = false;
    bool yLatAsc = false;
    bool yLatDesc = false;
    for(int j = 0; j < coordinates->height() - 1; j++) {
        for(int i = 0; i < coordinates->width(); i++) { // go through all x layers
            double2 p1 = coordinates->get(i, j);
            double2 p2 = coordinates->get(i, j + 1);
            if(p1.lon < p2.lon)yLonAsc = true;
            if(p1.lon > p2.lon)yLonDesc = true;
            if(p1.lat < p2.lat)yLatAsc = true;
            if(p1.lat > p2.lat)yLatDesc = true;
        }
    }
    
    // no orientation, if asc&desc values are both true for any lon/lat for x or y
    if((xLonAsc && xLonDesc) || (yLatAsc && yLatDesc))return OM_NONE;
    
    // primitives, created by rot90°
    if(xLonAsc && yLatAsc)return OM_NEMO;
    if(xLatAsc && yLonDesc)return OM_NEMO270CW;
    if(xLonDesc && yLatDesc)return OM_NEMO180CW;
    if(xLatDesc && yLonAsc)return OM_NEMO90CW;
    
    // primitives, after applying mirroring on x-axis
    if(xLonDesc && yLatAsc)return OM_NEMOMX;
    if(yLonDesc && xLatDesc)return OM_NEMO270CWMX;
    if(xLonAsc && yLatDesc)return OM_NEMO180CWMX;
    if(xLatAsc && yLonAsc)return OM_NEMO90CWMX;
    
    return -1;
}

int2 NCGrid::NEMOToOrientation(const int i, const int j, const int orientation) {
    // use here width/height functions that retrieve values in NEMO system of the loaded data!
    switch (orientation) {
        case OM_NEMO:
            // identity
            return int2(i, j);
            break;
        case OM_NEMOMX:
            return mirX(int2(i, j), width()); // mirror in x direction
            break;
        case OM_NEMO90CW:
            // rotate by 90° clockwise
            return rot90CW(int2(i, j), height()); // transformation function
            break;
        case OM_NEMO90CWMX:
        {
            // rotate first by 90° clockwise and mirror coords then
            int2 index = rot90CW(int2(i, j), height());
            return mirX(index, height()); // use here height, cause after 90° rotation width & height change!
            break;
        }
        case OM_NEMO180CW:
        {
            // rotate by 180° with specialized function(faster!)
            return rot180(int2(i, j), width(), height());
            break;
        }
        case OM_NEMO180CWMX:
        {
            int2 index = rot180(int2(i, j), width(), height());
            return mirX(index, width());
            break;
        }
        case OM_NEMO270CW:
        {
            // 270CW == 90CCW!
            return rot90CCW(int2(i, j), width());
            break;
        }
        case OM_NEMO270CWMX:
        {
            int2 index = rot90CCW(int2(i, j), width());
            // note that now height/width swapped!
            return mirX(index, height());
            break;
        }
        default:
        {
            std::cout<<">> error: unknown orientation"<<std::endl;
            return int2(-1, -1);
            break;
        }
    }
}

int2 NCGrid::OrientationToNEMO(const int i, const int j, const int orientation) {
    switch (orientation) {
        case OM_NEMO:
            // identity
            return int2(i, j);
            break;
        case OM_NEMOMX:
            return mirX(int2(i, j), coordinates->width()); // mirror in x direction
            break;
        case OM_NEMO90CW:
            // rotate by 90° counterclockwise
            return rot90CCW(int2(i, j), coordinates->width()); // transformation function
            break;
        case OM_NEMO90CWMX:
        {
            //mirror first in x-direction
            int2 index = mirX(int2(i, j), coordinates->width());
            return rot90CCW(index, coordinates->width()); // then rot 90° counter clockwise
            break;
        }
        case OM_NEMO180CW:
        {
            // rotate by 180° with specialized function(faster!)
            return rot180(int2(i, j), coordinates->width(), coordinates->height());
            break;
        }
        case OM_NEMO180CWMX:
        {
            //mirror first in x-direction
            int2 index = mirX(int2(i, j), coordinates->width());
            return rot180(index, coordinates->width(), coordinates->height()); // then rot 180°
            break;
        }
        case OM_NEMO270CW:
        {
            // 270CW == 90CCW!
            return rot90CW(int2(i, j), coordinates->height());
            break;
        }
        case OM_NEMO270CWMX:
        {
            //mirror first in x-direction
            int2 index = mirX(int2(i, j), coordinates->width());
            return rot90CW(index, coordinates->height()); // then 90° clockwise!
            break;
        }
        default:
        {
            std::cout<<">> error: unknown orientation"<<std::endl;
            return int2(-1, -1);
            break;
        }
    }
}

int NCGrid::width() {
    switch (orientation) {
        case OM_NEMO:
        case OM_NEMOMX:
        case OM_NEMO180CW:
        case OM_NEMO180CWMX:
            return coordinates->width();
            break;
        case OM_NEMO270CW:
        case OM_NEMO270CWMX:
        case OM_NEMO90CW:
        case OM_NEMO90CWMX:
            return coordinates->height();
        default:
            std::cout<<">> error: unknown orientation"<<std::endl;
            return -1;
            break;
    }
}

int NCGrid::height() {
    switch (orientation) {
        case OM_NEMO:
        case OM_NEMOMX:
        case OM_NEMO180CW:
        case OM_NEMO180CWMX:
            return coordinates->height();
            break;
        case OM_NEMO270CW:
        case OM_NEMO270CWMX:
        case OM_NEMO90CW:
        case OM_NEMO90CWMX:
            return coordinates->width();
        default:
            std::cout<<">> error: unknown orientation"<<std::endl;
            return -1;
            break;
    }
}

double2 NCGrid::getCoordinate(int i, int j) {
    // the indices are given in NEMO orientation --> reorient to the data's orientation
    int2 index = NEMOToOrientation(i, j, orientation);
    assert(coordinates);
    return coordinates->get(index.x, index.y);
}

int2 NCGrid::getOriginalIndices(int i, int j) {
    // the indices are given in NEMO orientation --> reorient to the data's orientation
    int2 index = NEMOToOrientation(i, j, orientation);
    return index;
}

double NCGrid::getDepth(int i, int j) {
    // the indices are given in NEMO orientation --> reorient to the data's orientation
    int2 index = NEMOToOrientation(i, j, orientation);
    assert(bathymetry);
    return bathymetry->get(index.x, index.y);
}

void NCGrid::setCoordinate(int i, int j, double2 val) {
    int2 index = NEMOToOrientation(i, j, orientation);
    assert(coordinates);
    coordinates->set(index.x, index.y, val);
}

void NCGrid::setDepth(int i, int j, double val) {
    int2 index = NEMOToOrientation(i, j, orientation);
    assert(bathymetry);
    bathymetry->set(index.x, index.y, val);
}
/// function to initiate coordinate & bathymetry array
void NCGrid::init(const int width, const int height) {
    SafeDelete(coordinates);
    SafeDelete(bathymetry);
    coordinates = new TArray<double2>(width, height);
    bathymetry = new TArray<double>(width, height);
    // set bathymetry to zero
    for(int i = 0; i < width * height; i++)bathymetry->set(i, 0);
}

void NCGrid::transpose() {
    // transpose arrays
    coordinates->transpose();
    bathymetry->transpose();
    

}