//WireCoilB.c

#include <math.h>

#if defined(WIN32)
#define SHARED_API __declspec(dllexport)
#else
#define SHARED_API
#endif

SHARED_API double dBx(int n, double* args){
    return (args[1] * sin(args[0]) + args[2] * cos(args[0]) + args[3]) / pow((args[4] + args[5] * cos(args[0]) + args[6] * sin(args[0])),1.5);}

SHARED_API double dBy(int n, double* args){
    return args[1] * cos(args[0]) / pow((args[2] + args[3] * cos(args[0]) + args[4] * sin(args[0])),1.5);}

SHARED_API double dBz(int n, double* args){
    return args[1] * sin(args[0]) / pow((args[2] + args[3] * cos(args[0]) + args[4] * sin(args[0])),1.5);}