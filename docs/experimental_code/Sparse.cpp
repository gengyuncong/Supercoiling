/*
 * author: max klein
 */
#include <iostream>
#include <vector>

#include "Sparse.h"

using std::cout;
using std::endl;
using std::vector;

typedef int Obs;
typedef vector<Obs> C;
typedef vector<C> CC;

#define ARRVEC(Elem, vec, arr) Elem vec(arr, arr + sizeof(arr) / sizeof(arr[0]))

int obsArr0[] = {1,0,1,3,2,2,0,4,3,1,1,2,4,4,1,1,1,3,2,2,0,4,0,4,1,4,4,2,3,3,2,2,0,2,0,3,2,1,4,2,0,0,3,3,4,4,2,0,0,0};
int obsArr1[] = {4,1,3,2,1,2,0,0,1,1,3,4,3,1,0,1,1,4,0,4,0,4,4,1,4,4,4,3,1,0,2,4,2,3,2,3,2,1,2,1,2,1,4,1,1,1,3,2,4,2};
int obsArr2[] = {2,0,3,1,1,3,1,4,3,3,2,3,1,4,0,0,3,4,0,3,4,4,2,2,4,3,1,1,0,4,0,2,1,4,2,0,0,1,1,2,2,3,4,3,1,0,1,3,1,4};
int obsArr3[] = {1,3,0,4,3,2,3,4,3,3,2,2,3,4,0,3,3,4,3,3,1,4,2,2,2,3,3,4,0,1,1,0,3,1,1,0,1,3,1,4,3,0,2,2,4,1,0,0,2,0};
ARRVEC(C, obsVec0, obsArr0);
ARRVEC(C, obsVec1, obsArr1);
ARRVEC(C, obsVec2, obsArr2);
ARRVEC(C, obsVec3, obsArr3);


int obsArr2D0[100][3] = {{0,2,0},{0,1,1},{0,2,1},{1,1,0},{1,2,1},{2,0,0},{0,0,1},{1,0,0},{1,1,1},{1,2,0},
                         {2,2,0},{0,1,0},{1,2,1},{1,1,0},{2,0,0},{0,2,2},{0,1,2},{1,1,0},{2,0,2},{1,2,1},
                         {0,0,1},{1,2,2},{0,0,1},{1,0,0},{0,2,1},{1,0,1},{2,1,2},{1,1,1},{0,2,0},{1,0,2},
                         {2,1,0},{2,1,1},{2,0,2},{0,1,2},{0,2,1},{0,1,2},{1,0,2},{2,2,0},{1,0,0},{2,0,0},
                         {1,1,1},{0,1,1},{1,2,2},{0,1,0},{1,1,2},{2,1,2},{0,2,0},{2,1,2},{0,0,1},{0,0,0},
                         {2,1,0},{0,2,0},{1,2,1},{1,0,0},{0,0,2},{2,0,1},{1,0,2},{2,2,1},{2,1,1},{1,2,2},
                         {2,2,1},{2,0,2},{0,0,1},{2,2,2},{0,0,2},{1,0,2},{2,2,0},{0,2,2},{2,1,0},{0,2,1},
                         {2,0,1},{1,1,2},{1,2,0},{2,1,2},{0,0,0},{0,2,0},{2,1,2},{1,2,0},{0,0,0},{0,0,0},
                         {2,1,0},{1,0,2},{0,2,1},{1,1,0},{0,1,2},{0,2,2},{0,1,1},{0,0,0},{2,1,1},{2,1,1},
                         {1,1,1},{0,2,2},{2,0,0},{2,1,1},{1,2,1},{0,1,0},{2,1,2},{0,2,1},{1,1,1},{0,0,1},};


namespace experimental {

void test_Sparse()
{
    cout << "test_Sparse" << endl;

    Sparse<Obs> sparse;
    sparse.addObs(obsVec0);

    cout << endl;
}

void test_Sparse_C()
{
    cout << "test_Sparse_C" << endl;

    CC obsVec2D0;
    for (int i=0; i<100; i++) {
        C obsVec(&obsArr2D0[i][0], &obsArr2D0[i][3]);
        obsVec2D0.push_back(obsVec);
    }

    Sparse<C> sparse2D;
    sparse2D.addObs(obsVec2D0);

//    cout << "Sparse.printObs()" << endl;
//    sparse2D.printObs();
//
//    cout << "Sparse.printVals()" << endl;
//    sparse2D.printVals();

    cout << "Sparse.printObsFlat()" << endl;
    sparse2D.printObsFlat();

    cout << endl;
}

}

int main()
{
//    experimental::test_Sparse();
    experimental::test_Sparse_C();
}
