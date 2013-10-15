//
//  main.cpp
//  gmetric
//
//  Created by Leonhard Spiegelberg on 17.09.13.
//  Copyright (c) 2013 Leonhard Spiegelberg. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <netcdfcpp.h>
#include <netcdf.h>
#include <vector>
#include <cmath>
#include <cstring>
#include "helper.h"
#include "Vector2.h"
#include "Vector3.h"
#include "Vector2Test.h"
#include "Vector3Test.h"
#include "helperTest.h"
#include "MetricsTest.h"
#include "NCGrid.h"
#include "NCGridTest.h"
#include "AreaCalculator.h"
#include "AreaCalculatorTest.h"
#include "VolumeCalculator.h"
#include "VolumeCalculatorTest.h"

#define INTEGRATIONMETHOD_LINEAR 0x10
#define INTEGRATIONMETHOD_SPLINE 0x20
#define INTEGRATIONMODE_NONE 0x0
#define INTEGRATIONMODE_AREA 0x1
#define INTEGRATIONMODE_VOLUME 0x2

// function prototypes

/// little helper to check coords in degrees
bool checkIfCoordsInRange(const double2& P) {
    using namespace std;
    bool error = false;
    if(P.x < -180.0)error = true;
    if(P.x > 180.0)error = true;
    if(P.y < -90.0)error = true;
    if(P.y > 90.0)error = true;
    
    if(error) {
        cout<<">> error: invalid coordinate found: "<<P<<endl;
    }
    return !error;
}

/// prints short and precise how to use the program
void printShortUsage() {
    using namespace std;
    cout<<"usage: gmetric [--help] [--test] [--area lon1 lat1 lon2 lat2 | --area <file> | --volume lon1 lat1 lon2 lat2 lon3 lat3] [--output-segments <file>] [--output-areas <file>] [--bathymetry-var varname] [--lon-var varname] [--lat-var varname] [--linear | --spline] [--realworld] input_file"<<endl;
    cout<<endl;
    cout<<"\tfor help call gmetric --help"<<endl;
}

/// prints long form of usage
void printUsage() {
    using namespace std;
    cout<<"(c) 2013 L.Spiegelberg"<<endl;
    printShortUsage();
    cout<<endl;
    cout<<"list of commands:"<<endl;
    cout<<"\t--help\t\t\t\t\tprints help"<<endl;
    cout<<"\t--test\t\t\t\t\trun test cases"<<endl;
    cout<<"\t--area lon1 lat1 lon2 lat2 \tintegrate over line with\n\t\t\t\t\t\t\tcoordinates (lon1, lat1) and\n\t\t\t\t\t\t\t(lon2, lat2) in degrees"<<endl;
    cout<<"\t--area <file>\t\t\t\tintegrate over line with input \n\t\t\t\t\t\t\tfile containing line segments\n\t\t\t\t\t\t\t\n\t\t\t\t\t\t\tsyntax of input file has to be\n\t\t\t\t\t\t\t# this is a comment line\n\t\t\t\t\t\t\t# lon1 lat1 lon2 lat2 (in degrees)\n\t\t\t\t\t\t\t   60.3 10.4 20.0 10.0\n\t\t\t\t\t\t\t  .... .... .... ....\n\t\t\t\t\t\t\t  ...."<<endl<<endl;
    cout<<"\t--volume lon1 lat1\t\t"<<"\tintegrate over triangle defined\n\t         lon2 lat2 \t\tby coordinates in degrees\n\t         lon3 lat3"<<endl;
    cout<<"\tinput_file\t\t\t\tnetCDF file of the input grid"<<endl;
    cout<<"\t--output-segments <file>\toutput line segments to file"<<endl;
    cout<<"\t--output-areas    <file>\toutput covered areas to file"<<endl;
    cout<<"\t--bathymetry-var varname\tspecify how netCDF variable of \n\t\t\t\t\t\t\tbathymetry is called in the \n\t\t\t\t\t\t\tinput file\n\t\t\t\t\t\t\t(default is bathymetry)"<<endl;
    cout<<"\t--lon-var varname\t\t\tspecify how netCDF variable of \n\t\t\t\t\t\t\tlongitude is called in the \n\t\t\t\t\t\t\tinput file(default is lon)"<<endl;
    cout<<"\t--lat-var varname\t\t\tspecify how netCDF variable of \n\t\t\t\t\t\t\tlatitude is called in the \n\t\t\t\t\t\t\tinput file(default is lat)"<<endl;
    cout<<"\t--realworld\t\t\t\tcompute real world data where \n\t\t\t\t\t\t\tthe shape of the earth is based \n\t\t\t\t\t\t\ton a sphere with\n\t\t\t\t\t\t\tR := 6378.137km"<<endl;
    cout<<"\t--linear\t\t\t\t\tcell's edges are approximated \n\t\t\t\t\t\t\twith lines"<<endl;
    cout<<"\t--spline\t\t\t\t\tcell's edges are approximated \n\t\t\t\t\t\t\twith splines"<<endl;
}


// CPPUNIT defines for tests
//CPPUNIT's ugly static variable hiding concept...
CPPUNIT_TEST_SUITE_REGISTRATION(VolumeCalculatorTest);
/*CPPUNIT_TEST_SUITE_REGISTRATION(AreaCalculatorTest);
CPPUNIT_TEST_SUITE_REGISTRATION(NCGridTest);
CPPUNIT_TEST_SUITE_REGISTRATION(Vector2Test);
CPPUNIT_TEST_SUITE_REGISTRATION(Vector3Test);
CPPUNIT_TEST_SUITE_REGISTRATION(helperTest);*/

// helper to run tests
bool runTests() {
    using namespace std;
    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
    runner.addTest( registry.makeTest() );
    bool wasSuccessful = runner.run( "", false );
    
    if(wasSuccessful) {
        cout<<" >> all tests succeeded! "<<endl;
    }
    else cout<<" >> test suite failed! "<<endl;
    return wasSuccessful;
}

// output will be written to vector v
// only linear supported yet!
double subdivide_volume(NCGrid *grid, std::string areasFile, std::vector<double>& v, const double2& A, const double2& B, const double2& C, const bool realWorldData, int level) {
    // if level <= 1, return direct volume
    if(level <= 1) {
        double volume = VolumeCalculator::instance().linear(grid, A, B, C, areasFile, realWorldData, EARTH_RADIUS);
        
        // output something
        std::cout<<"part result "<<volume<<" km^3..."<<std::endl;
        
        // push pack to vector
        v.push_back(volume);
        
        return volume;
    }
    else {
        double volume = 0;
        // subdivide!
        double2 ACM = (A + C) / 2.0;
        double2 ABM = (A + B) / 2.0;
        double2 BCM = (B + C) / 2.0;
        
        
        // subdivide in 4 triangles!
        volume += subdivide_volume(grid, areasFile, v, A, ACM, ABM, realWorldData, level - 1);
        volume += subdivide_volume(grid, areasFile, v, ACM, ABM, BCM, realWorldData, level - 1);
        volume += subdivide_volume(grid, areasFile, v, ABM, BCM, B, realWorldData, level - 1);
        volume += subdivide_volume(grid, areasFile, v, C, ACM, BCM, realWorldData, level - 1);
        return volume;
    }
}

int main(int argc, const char * argv[])
{
    using namespace std;
    // all variables, set to default first
    int     integrationMethod       = INTEGRATIONMETHOD_LINEAR;
    int     integrationMode         = INTEGRATIONMODE_NONE;
    string  bathymetryVarName    = "bathymetry";
    string  lonVarName           = "lon";
    string  latVarName           = "lat";
    bool    computeRealWorldData   = false;
    bool    outputSegments         = false;
    bool    outputAreas            = false;
    string  segmentFile = "";
    string  areasFile;
    string  ncFilePath = "none";
    string  coordListFileArea = "";     // use instead of plain coords a file?
    string  coordListFileVolume = "";   // use instead of plain coords a file?
    double2 lA, lB; // line coords
    double2 tA, tB, tC; // triangle coords
    
    bool swapDims = false;          // undocument feature to handle wrong ncFiles
    string  partsFile = "";         // if a file is given, write data of it to output
    int subdivide_level = 0;        // level for which triangle should be subdivided
    string subdivide_file = "";     // file where to save volumes of subdivides!
    
    // if no arguments is given, print usage
    if (argc == 1) {
        printShortUsage();
        return 0;
    }
    
    // collect parameters
    for(int i = 1; i < argc;i++) {
        if(0 == strcmp(argv[i], "--help")) {
            printUsage();
            return 0;
        }
        else
        if(0 == strcmp(argv[i], "--test")) {
            runTests();
            return 0;
        }
        else
        if(0 == strcmp(argv[i], "--bathymetry-var")) {
            // assert there are enough arguments left
            if(i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --bathymetry-var"<<endl;
                return 0;
            }
            i++;
            bathymetryVarName = argv[i];
        }
        else
        if(0 == strcmp(argv[i], "--lon-var")) {
            // assert there are enough arguments left
            if(i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --lon-var"<<endl;
                return 0;
            }
            i++;
            lonVarName = argv[i];
        }
        else
        if(0 == strcmp(argv[i], "--lat-var")) {
            // assert there are enough arguments left
            if(i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --lat-var"<<endl;
                return 0;
            }
            i++;
            latVarName = argv[i];
        }
        else
        if(0 == strcmp(argv[i], "--realworld"))computeRealWorldData = true;
        else if(0 == strcmp(argv[i], "--linear"))integrationMethod = INTEGRATIONMETHOD_LINEAR;
        else if(0 == strcmp(argv[i], "--spline"))integrationMethod = INTEGRATIONMETHOD_SPLINE;
        else if(0 == strcmp(argv[i], "--output-segments")) {
            // assert there are enough arguments left
            if(i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --output-segments"<<endl;
                return 0;
            }
            i++;
            outputSegments = true;
            segmentFile = argv[i];
        }
        else if(0 == strcmp(argv[i], "--output-areas")) {
            // assert there are enough arguments left
            if(i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --output-areas"<<endl;
                return 0;
            }
            i++;
            outputAreas = true;
            areasFile = argv[i];
        }
        else if(0 == strcmp(argv[i], "--parts")) {
            // assert there are enough arguments left
            if(i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --parts"<<endl;
                return 0;
            }
            i++;
            partsFile = argv[i];
        }
        // read input vars
        else if(0 == strcmp(argv[i], "--area")) {
            integrationMode |= INTEGRATIONMODE_AREA;
            // assert there are enough arguments left
            if(i + 4 >= argc || i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --area"<<endl;
                return 0;
            }
            
            // decide on next argument if a file is loaded or coords are used!
            bool loadfile = false;
            if(i + 2 == argc)loadfile = true;
            if(i + 2 < argc)if(argv[i + 2][0] == '-')loadfile = true;
            
            if(loadfile) {
                coordListFileArea = argv[++i];
            } else {
                lA.x = atof(argv[++i]);
                lA.y = atof(argv[++i]);
                lB.x = atof(argv[++i]);
                lB.y = atof(argv[++i]);
            }
        }
        else if(0 == strcmp(argv[i], "--volume")) {
            integrationMode |= INTEGRATIONMODE_VOLUME;
            
            // assert there are enough arguments left
            if(i + 6 >= argc || i + 1 >= argc) {
                cout<<">> error: invalid number of arguments for --volume"<<endl;
                return 0;
            }
            
            // decide on next argument if a file is loaded or coords are used!
            bool loadfile = false;
            if(i + 2 == argc)loadfile = true;
            if(i + 2 < argc)if(argv[i + 2][0] == '-')loadfile = true;
            
            if(loadfile) {
                coordListFileVolume = argv[++i];
            } else {
                tA.x = atof(argv[++i]);
                tA.y = atof(argv[++i]);
                tB.x = atof(argv[++i]);
                tB.y = atof(argv[++i]);
                tC.x = atof(argv[++i]);
                tC.y = atof(argv[++i]);
            }
        }
        else if(0 == strcmp(argv[i], "--subdivide")) {
            // assert there are enough arguments left
            if(i + 2 >= argc) {
                cout<<">> error: invalid number of arguments for --subdivide"<<endl;
                return 0;
            }
            
            subdivide_level = atoi(argv[++i]);
            subdivide_file = argv[++i];
            assert(subdivide_level > 0);
        }
        else if(0 == strcmp(argv[i], "--swapdims"))
        {
            swapDims = true;
        }
        else {
            // this has to be the input file...
            ncFilePath = argv[i];
        }
    }
    
    // no integration mode defined?
    if(INTEGRATIONMODE_NONE == integrationMode) {
        cout<<">> error: you must define at least one integration mode(--area ... or --volume ...)"<<endl;
        return 0;
    }
    
    // no input file?
    if(strcmp(ncFilePath.c_str(), "none") == 0) {
        cout<<">> error: you must specify a netCDF input file"<<endl;
        return 0;
    }
    
    // now read first the input file in...
    NCGrid *ncGrid = new NCGrid;
    if(ncGrid->open(ncFilePath, bathymetryVarName, lonVarName, latVarName, swapDims) < 0) {
        cout<<">> error: netCDF file could not be opened sucessfully"<<endl;
        SafeDelete(ncGrid);
        return 0;
    }
    
    // set cout precision to 4 digits
    cout.precision(4);
    cout<<fixed;
    
    
    // now calculate depending on selected mode metric
    if(INTEGRATIONMODE_AREA == integrationMode) {
        double area = 0;
        cout<<">> info: starting computation..."<<endl;
        // use a list?
        if(coordListFileArea.length() == 0) {
            if(!checkIfCoordsInRange(lA))return 0;
            if(!checkIfCoordsInRange(lB))return 0;
            
            if(INTEGRATIONMETHOD_LINEAR == integrationMethod)
                area = AreaCalculator::instance().linear(ncGrid, lA, lB, segmentFile, computeRealWorldData, EARTH_RADIUS);
            else if(INTEGRATIONMETHOD_SPLINE == integrationMethod)
                area = AreaCalculator::instance().spline(ncGrid, lA, lB, segmentFile, computeRealWorldData, EARTH_RADIUS);
        }
        else {
            
            // if parts file is valid, remove it
            if(partsFile.length() != 0) {
                remove(partsFile.c_str());
                ofstream ofs(partsFile.c_str());
                ofs<<"# lon1 lat1 lon2 lat2 distance area"<<endl;
                ofs.close();
            }
            
            // if output segments is wished, clear existing file!
            if(outputSegments) {
                remove(segmentFile.c_str());
            }
            
            // open file and read coords, sum up area!
            ifstream ifs(coordListFileArea.c_str(), ios::in);
            char buffer[1024];
            // errors?
            if(!ifs) {
                cout<<">> error: failed opening file "<<coordListFileArea<<endl;
                return 0;
            }
            
            int line = 0;
            // easy format, scan lines and retrieve lon1 lat1 lon2 lat2
            while(ifs.getline(buffer, 1024)) {
                
                line++;
                
                // first character # sign? => skip!
                if(buffer[0] == '#')continue;
                
                // split string in token
                char *ptr = NULL;
                char delimiter[] = ",; ";
                // initialize
                ptr = strtok(buffer, delimiter);
                
                int pos = 0;
                double v[4];
                while(ptr != NULL && pos < 4) {
                    // current token it ptr
                    v[pos++] = atof(ptr);
                    
                    // get next token
                    ptr = strtok(NULL, delimiter);
                }
                
                // set points
                lA = double2(v[0], v[1]);
                lB = double2(v[2], v[3]);
                
                if(!checkIfCoordsInRange(lA))return 0;
                if(!checkIfCoordsInRange(lB))return 0;
                
                if(pos != 4) {
                    cout<<">> error: not enough input values in line "<<line<<endl;
                    return 0;
                } else {
                    
                    // calc metric add to area
                    double partarea = 0;
                    if(INTEGRATIONMETHOD_LINEAR == integrationMethod)
                        partarea = AreaCalculator::instance().linear(ncGrid, lA, lB, segmentFile, computeRealWorldData, EARTH_RADIUS);
                    else if(INTEGRATIONMETHOD_SPLINE == integrationMethod)
                        partarea = AreaCalculator::instance().spline(ncGrid, lA, lB, segmentFile, computeRealWorldData, EARTH_RADIUS);
                    
                    // if partsfile is valid write data to it...
                    if(partsFile.length() != 0) {
                        ofstream ofs(partsFile.c_str(), ios::app);
                        double partdistance = AreaCalculator::instance().llDistanceUnitSphere(lA, lB, computeRealWorldData) * EARTH_RADIUS;
                        ofs<<lA.x<<" "<<lA.y<<" "<<lB.x<<" "<<lB.y<<" "<<partdistance<<" "<<partarea<<endl;
                    }
                    
                    
                    // print out current result
                    cout<<">> info: part result "<<partarea<<"km^2"<<endl;
                    area += partarea;
                }
            }
            
            // errors during reading?
            if(ifs.bad()) {
                cout<<">> error: failure during reading "<<coordListFileArea<<endl;
                return 0;
            }
            // close stream
            ifs.close();
        }
        
        // print out result
        cout<<">> result: "<<area<<" km^2"<<endl;
    }


    
    // now calculate depending on selected mode metric (M_2 <-> volume)
    if(INTEGRATIONMODE_VOLUME == integrationMode) {
        double volume = 0;
        cout<<">> info: starting computation..."<<endl;
        // use subdivision?
        if(subdivide_level != 0 && subdivide_file.length() != 0) {
            // vector to save data
            vector<double>  volumes;
            ofstream ofs(subdivide_file.c_str());
            // errors?
            if(!ofs) {
                cout<<">> error: failed opening file "<<subdivide_file<<endl;
                return 0;
            }
            
            // now calc all subdivisions
            volume = subdivide_volume(ncGrid, areasFile, volumes, tA, tB, tC, computeRealWorldData, subdivide_level);
            
            if(!volumes.empty()) {
                ofs<<"# volume of subdivided triangle"<<endl;
                ofs<<fixed;
                ofs.precision(8);
                for(vector<double>::const_iterator it = volumes.begin(); it != volumes.end(); ++it) {
                    ofs<<*it<<endl;
                }
            }
            
            ofs.close();
            
            cout<<"computation finished, volume = "<<volume<<" km^3"<<endl;
            volume = VolumeCalculator::instance().linear(ncGrid, tA, tB, tC, areasFile, computeRealWorldData, EARTH_RADIUS);
            cout<<"computation finished, volume = "<<volume<<" km^3"<<endl;
        }
        
        // use a list?
        else if(coordListFileVolume.length() == 0) {
            if(!checkIfCoordsInRange(tA))return 0;
            if(!checkIfCoordsInRange(tB))return 0;
            if(!checkIfCoordsInRange(tC))return 0;
            
            if(INTEGRATIONMETHOD_LINEAR == integrationMethod)
                volume = VolumeCalculator::instance().linear(ncGrid, tA, tB, tC, areasFile, computeRealWorldData, EARTH_RADIUS);
            else if(INTEGRATIONMETHOD_SPLINE == integrationMethod)
                volume = VolumeCalculator::instance().spline(ncGrid, tA, tB, tC, areasFile, computeRealWorldData, EARTH_RADIUS);
        }
        else {
            // if output segments is wished, clear existing file!
            if(outputAreas) {
                remove(areasFile.c_str());
            }
            
            // open file and read coords, sum up area!
            ifstream ifs(coordListFileVolume.c_str(), ios::in);
            char buffer[1024];
            // errors?
            if(!ifs) {
                cout<<">> error: failed opening file "<<coordListFileVolume<<endl;
                return 0;
            }
            
            int line = 0;
            // easy format, scan lines and retrieve lon1 lat1 lon2 lat2
            while(ifs.getline(buffer, 1024)) {
                
                line++;
                
                // first character # sign? => skip!
                if(buffer[0] == '#')continue;
                
                // split string in token
                char *ptr = NULL;
                char delimiter[] = ",; ";
                // initialize
                ptr = strtok(buffer, delimiter);
                
                int pos = 0;
                double v[6];
                while(ptr != NULL && pos < 6) {
                    // current token it ptr
                    v[pos++] = atof(ptr);
                    
                    // get next token
                    ptr = strtok(NULL, delimiter);
                }
                
                // set points
                tA = double2(v[0], v[1]);
                tB = double2(v[2], v[3]);
                tC = double2(v[2], v[3]);
                
                if(!checkIfCoordsInRange(tA))return 0;
                if(!checkIfCoordsInRange(tB))return 0;
                if(!checkIfCoordsInRange(tC))return 0;
                
                if(pos != 6) {
                    cout<<">> error: not enough input values in line "<<line<<endl;
                    return 0;
                } else {
                    
                    // calc metric add to area
                    double partvolume = 0;
                    if(INTEGRATIONMETHOD_LINEAR == integrationMethod)
                        partvolume = volume = VolumeCalculator::instance().linear(ncGrid, tA, tB, tC, areasFile, computeRealWorldData, EARTH_RADIUS);
                    else if(INTEGRATIONMETHOD_SPLINE == integrationMethod)
                        partvolume = volume = VolumeCalculator::instance().spline(ncGrid, tA, tB, tC, areasFile, computeRealWorldData, EARTH_RADIUS);
                    
                    // print out current result
                    cout<<">> info: part result "<<partvolume<<"km^3"<<endl;
                    volume += partvolume;
                }
            }
            
            // errors during reading?
            if(ifs.bad()) {
                cout<<">> error: failure during reading "<<coordListFileVolume<<endl;
                return 0;
            }
            // close stream
            ifs.close();
        }
        
        // print out result
        cout<<">> result: "<<volume<<" km^3"<<endl;
    }
    
    
    // clean mem
    SafeDelete(ncGrid);
    return 0;
}

