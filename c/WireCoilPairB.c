//WireCoilPairB.c

#include <math.h>
using namespace std;

double lfdBx(int n, double args[n]){
    return args[1] * (-args[2]*args[5]*sin(args[0]) - args[2]*args[4]*cos(args[0]) + pow(args[2],2)) / pow((pow(args[3] + args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);
}

double rtdBx(int n, double args[n]){
    return args[1] * (-args[2]*args[5]*sin(args[0]) - args[2]*args[4]*cos(args[0]) + pow(args[2],2)) / pow((pow(args[3] - args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);
}

//double cst, double R, double px, double py, double pz, double d
// args[1]      args[2]   args[3]   args[4]     args[5]   args[6]

double lfdBy(int n, double args[n]){
    return args[1] * (args[2]*(args[3]+args[6])*cos(args[0])) / pow((pow(args[3] + args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);
}

double rtdBy(int n, double args[n]){
    return args[1] * (args[2]*(args[3]-args[6])*cos(args[0])) / pow((pow(args[3] - args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);
}

double lfdBz(int n, double args[n]){
    return args[1] * (args[2]*(args[3]+args[6])*sin(args[0])) / pow((pow(args[3] + args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);
}
double rtdBz(int n, double args[n]){
    return args[1] * (args[2]*(args[3]-args[6])*sin(args[0])) / pow((pow(args[3] - args[6],2) + pow(args[4]-args[2]*cos(args[0]),2) + pow(args[5]-args[2]*sin(args[0]),2)),1.5);
}

//lfdBx = lambda a: self.cst*(-self.R*ppr[2]*sin(a) - self.R*ppr[1]*cos(a) + self.R**2) / (((ppr[0] + self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - self.R*sin(a))**2)**(3/2))

//rtdBx = lambda a: self.cst*(-self.R*ppr[2]*sin(a) - self.R*ppr[1]*cos(a) + self.R**2) / (((ppr[0] - self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - self.R*sin(a))**2)**(3/2))

//lfdBy = lambda a: self.cst*(self.R*cos(a)*(ppr[0] + self.d)) / (((ppr[0] + self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - self.R*sin(a))**2)**(3/2))

//rtdBy = lambda a: self.cst*(self.R*cos(a)*(ppr[0] - self.d)) / (((ppr[0] - self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - self.R*sin(a))**2)**(3/2))

//lfdBz = lambda a: self.cst*(self.R*sin(a)*(ppr[0] + self.d)) / (((ppr[0] + self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - self.R*sin(a))**2)**(3/2))

//rtdBz = lambda a: self.cst*(self.R*sin(a)*(ppr[0] - self.d)) / (((ppr[0] - self.d)**2 + (ppr[1] - self.R*cos(a))**2 + (ppr[2] - self.R*sin(a))**2)**(3/2))