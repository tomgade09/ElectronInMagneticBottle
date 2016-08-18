//WireCoilPairB.c

#include <math.h>
//using namespace std;

double lfdBx(int n, double args[n]){
    return args[1] * (-args[2]*args[5]*sin(args[0]) - args[2]*args[4]*cos(args[0]) + pow(args[2],2)) / pow((pow(args[3] + args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);}

double rtdBx(int n, double args[n]){
    return args[1] * (-args[2]*args[5]*sin(args[0]) - args[2]*args[4]*cos(args[0]) + pow(args[2],2)) / pow((pow(args[3] - args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);}

double lfdBy(int n, double args[n]){
    return args[1] * (args[2]*(args[3]+args[6])*cos(args[0])) / pow((pow(args[3] + args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);}

double rtdBy(int n, double args[n]){
    return args[1] * (args[2]*(args[3]-args[6])*cos(args[0])) / pow((pow(args[3] - args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);}

double lfdBz(int n, double args[n]){
    return args[1] * (args[2]*(args[3]+args[6])*sin(args[0])) / pow((pow(args[3] + args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);}

double rtdBz(int n, double args[n]){
    return args[1] * (args[2]*(args[3]-args[6])*sin(args[0])) / pow((pow(args[3] - args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);}