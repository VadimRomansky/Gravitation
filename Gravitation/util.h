#ifndef UTIL_H
#define UTIL_H

#include <string>

const int randomParameter = 1024;
const double pi = 4*atan(1.0);

double uniformDistribution();
double sqr(const double& a);
double cube(const double& a);
double min(const double& a, const double& b);
double min3(const double& a, const double& b, const double& c);
double min4(const double& a, const double& b, const double& c, const double& d);
double max(const double& a, const double& b);
double evaluateAngle(const double& ax, const double& ay, const double& az, const double& bx, const double& by, const double& bz);
std::string convertIntToString(int a);

#endif