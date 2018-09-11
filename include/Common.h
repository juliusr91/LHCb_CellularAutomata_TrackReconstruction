//
// Created by julius on 3/5/18.
//

#include <glob.h>
#include <vector>
#include <stdint-gcc.h>
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include <string>
#include <fstream>
#include <chrono>
#include <array>

#define NUMBER_OF_SENSORS 52
#define NEIGHBOUR_INDEX 2

#define MAX_SLOPE_X 0.6f
#define MAX_SLOPE_Y 0.6f

#define MIN_LENGTH_TRACK 3
#define MAX_PERCENTAGE_OF_USED_HITS 0.26
#define PRINT_BINARY true
#define MAX_RAD_DIFF 0.125f

#define NUMBER_LAYERS_DOUBLETS 50
#define MAXNUM_NEIGHBOURS_POINT 10 //for finding doublets
#define MAXNUM_NEIGHBOURS_DOUBLET 3 //for making the actual neighbours of doublets
#define MAXNUM_TRACKS 1000


//original event information object
struct EventInfo {
    size_t size;
    uint32_t numberOfModules;
    uint32_t numberOfHits;
    float* module_Zs;
    uint32_t* module_hitStarts;
    uint32_t* module_hitNums;
    uint32_t* hit_IDs;
    float* hit_Xs;
    float* hit_Ys;
    float* hit_Zs;
    std::vector<float> phi;


    EventInfo() = default;
    explicit EventInfo(const std::vector<uint8_t> &event){
        uint8_t * input = (uint8_t*) event.data();

        numberOfModules  = *((uint32_t*)input); input += sizeof(uint32_t);
        numberOfHits     = *((uint32_t*)input); input += sizeof(uint32_t);
        module_Zs        = (float*)input; input += sizeof(float) * numberOfModules;
        module_hitStarts = (uint32_t*)input; input += sizeof(uint32_t) * numberOfModules;
        module_hitNums   = (uint32_t*)input; input += sizeof(uint32_t) * numberOfModules;
        hit_IDs          = (uint32_t*)input; input += sizeof(uint32_t) * numberOfHits;
        hit_Xs           = (float*)  input; input += sizeof(float)   * numberOfHits;
        hit_Ys           = (float*)  input; input += sizeof(float)   * numberOfHits;
        hit_Zs           = (float*)  input; input += sizeof(float)   * numberOfHits;

        size = input - (uint8_t*) event.data();

        phi.reserve(numberOfHits);
        for (int eventNum = 0; eventNum < numberOfHits; eventNum++){
            phi.push_back(std::atan2(hit_Ys[eventNum], hit_Xs[eventNum]));
        }
    }
};

//neighbour object.
//can probably get rid of unless you want to include the skipped sensors.
struct Neighbour{
    int module = -1;
    int doubletNum = -1;

    Neighbour() = default;
    explicit Neighbour(int input1, int input2){
        module = input1;
        doubletNum = input2;
    }
};


//track object
struct Track{
    int length = 0;
    std::vector<int> listHits;
    float chi2 = 0;

    Track()= default;
    Track(int hit1, int hit2){
        listHits.reserve(5); //good parameter??
        listHits.emplace_back(hit1);
        listHits.emplace_back(hit2);
    }

    //copy instructor
    Track(const Track &obj){
        listHits = obj.listHits;
        length = obj.length;
        chi2 = obj.chi2;
    }
};

//the doublet struct.
//could probably be simplified by not taking into account phi?
struct Doublet{
    int hit1;
    int hit2;
    float phi1;
    float phi2;

    int numNeighbours = 0;
    std::array<Neighbour, MAXNUM_NEIGHBOURS_DOUBLET> neighbours;

    uint16_t status {1};
    uint16_t newStatus {1};

    Doublet()=default;
    explicit Doublet(int input1, int input2, float input3, float input4){
        hit1 = input1;
        hit2 = input2;
        phi1 = input3;
        phi2 = input4;
        numNeighbours = 0;
    }
    void makeNewNeighbour(Neighbour newNeighbour){
        neighbours[numNeighbours] = newNeighbour;
        numNeighbours++;
    }


};

//sorting by length and chi2
struct sortingOnLengthAndChi2 {
    inline bool operator() (const Track &track1, const Track &track2){
        if (track1.length != track2.length){
            return track1.length > track2.length;
        }
        return track1.chi2 < track2.chi2;

    }
};



//used in the lower bound search of the neighbour finding function
//could probably use the same for searching in the make doublets phase.
class search_comparator {
public:
    bool operator()(const Doublet &a, const int &b) const
    {
        return a.hit2 < b;
    }

    bool operator()(const int &a, const Doublet &b) const
    {
        return a < b.hit2;
    }
};




//the original code for the vectorised area calculation by Gerhard Raven can be found here:
//https://godbolt.org/g/bu2bfK
//the mathematical explanations can be found here:
//https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle

template <typename T,unsigned N> struct vec {
#ifdef __clang__
    typedef T type __attribute__ ((ext_vector_type( N  ))); // vector of N Ts
#else
    typedef T type __attribute__ ((vector_size( N * sizeof(T) ))); // vector of N Ts
#endif
};


template <typename T> struct Point
{
    T x,y,z;
    friend  Point operator-(const Point& a, const Point& b)
    { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
    friend  Point cross(const Point& a , const Point& b)
    { return { a.y*b.z-a.z*b.y, -a.x*b.z+a.z*b.x, a.x*b.y-a.y*b.x}; }
    T mag2() const { return x*x+y*y+z*z;  }
};

template <typename T> T two_area_sq_v( const Point<T>& a, const Point<T>& b, const Point<T>& c) {
    return cross(b-a,c-a).mag2() ;
}

using v8 = vec<float,8>::type;
template v8 two_area_sq_v(const Point<v8>& ,const Point<v8>& , const Point<v8>&);


template <typename T>
T two_area_sq(  Point<T> a,  Point<T> b,  Point<T> c) {
    using v4 = typename vec<T,4>::type;
    const auto va = v4{ a.x, a.y, a.z, 0 };
    const auto ba = v4{ b.x, b.y, b.z, 0 }-va;
    const auto ca = v4{ c.x, c.y, c.z, 0 }-va;

    const v4 x[] = { ba*ca[0], ba*ca[1], ba*ca[2] };

    const auto i = x[2][1]-x[1][2];
    const auto j = x[2][0]-x[0][2];
    const auto k = x[1][0]-x[0][1];

    return i*i+j*j+k*k;
}
template float two_area_sq( Point<float>,  Point<float>, Point<float>);
