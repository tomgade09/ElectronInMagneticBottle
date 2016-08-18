//WireCoilPairB.c

#include <math.h>
//using namespace std;

double lfdBx(int n, double args[n]){
    return args[1] * (-args[2]*args[5]*math::sin(args[0]) - args[2]*args[4]*math::cos(args[0]) + math::pow(args[2],2)) / math::pow((math::pow(args[3] + args[6],2) + math::pow(args[4]-args[2]*math::cos(args[0]),2) + math::pow(args[5]-args[2]*math::sin(args[0]),2)),1.5);}

double rtdBx(int n, double args[n]){
    return args[1] * (-args[2]*args[5]*math::sin(args[0]) - args[2]*args[4]*math::cos(args[0]) + math::pow(args[2],2)) / math::pow((math::pow(args[3] - args[6],2) + math::pow(args[4]-args[2]*math::cos(args[0]),2) + math::pow(args[5]-args[2]*math::sin(args[0]),2)),1.5);}

double lfdBy(int n, double args[n]){
    return args[1] * (args[2]*(args[3]+args[6])*math::cos(args[0])) / math::pow((math::pow(args[3] + args[6],2) + math::pow(args[4]-args[2]*math::cos(args[0]),2) + math::pow(args[5]-args[2]*math::sin(args[0]),2)),1.5);}

double rtdBy(int n, double args[n]){
    return args[1] * (args[2]*(args[3]-args[6])*math::cos(args[0])) / math::pow((math::pow(args[3] - args[6],2) + math::pow(args[4]-args[2]*math::cos(args[0]),2) + math::pow(args[5]-args[2]*math::sin(args[0]),2)),1.5);}

double lfdBz(int n, double args[n]){
    return args[1] * (args[2]*(args[3]+args[6])*math::sin(args[0])) / math::pow((math::pow(args[3] + args[6],2) + math::pow(args[4]-args[2]*math::cos(args[0]),2) + math::pow(args[5]-args[2]*math::sin(args[0]),2)),1.5);}

double rtdBz(int n, double args[n]){
    return args[1] * (args[2]*(args[3]-args[6])*math::sin(args[0])) / math::pow((math::pow(args[3] - args[6],2) + math::pow(args[4]-args[2]*math::cos(args[0]),2) + math::pow(args[5]-args[2]*math::sin(args[0]),2)),1.5);}