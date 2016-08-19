//WireCoilPairB.c

#include <math.h>

double dBx(int n, double args[n]){
    return (args[1] * sin(args[0]) + args[2] * cos(args[0]) + args[3]) / pow((args[4] + args[5] * cos(args[0]) + args[6] * sin(args[0])),1.5);}

double dBy(int n, double args[n]){
    return args[1] * cos(args[0]) / pow((args[2] + args[3] * cos(args[0]) + args[4] * sin(args[0])),1.5);}

double dBz(int n, double args[n]){
    return args[1] * sin(args[0]) / pow((args[2] + args[3] * cos(args[0]) + args[4] * sin(args[0])),1.5);}