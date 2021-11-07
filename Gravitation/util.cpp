#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <omp.h>


#include "util.h"

double uniformDistribution() {
	return (rand() % randomParameter + 0.5) / randomParameter;
}

double sqr(const double& a){
	return a*a;
}

double cube(const double& a){
	return a*a*a;
}

double min(const double& a, const double& b){
	if(a < b) {
		return a;
	} else {
		return b;
	}
}

double min3(const double& a, const double& b, const double& c){
	if(a < b) {
		return min(a, c);
	} else {
		return min(b, c);
	}
}

double min4(const double& a, const double& b, const double& c, const double& d){
	if(a < b) {
		if(c < d){
			return min(a, c);
		} else {
			return min(a, d);
		}
	} else {
		if(c < d){
			return min(b, c);
		} else {
			return min(b, d);
		}
	}
}

double evaluateAngle(const double& ax, const double& ay, const double& az, const double& bx, const double& by, const double& bz){
	double a = sqrt(ax*ax + ay*ay + az*az);
	double b = sqrt(bx*bx + by*by + bz*bz);

	double scalar = ax*bx + ay*by + az*bz;
	double cosab = scalar/(a*b);

	double cx = ay*bz - az*by;
	double cy = az*bx - ax*bz;
	double cz = ax*by - ay*bx;

	if(cz > 0){
		return 180*acos(cosab)/pi;
	} else {
		return -180*acos(cosab)/pi;
	}

	/*double c = sqrt(cx*cx + cy*cy + cz*cz);
	double sinab = c/(a*b);*/
}

double max(const double& a, const double& b){
	if(a > b) {
		return a;
	} else {
		return b;
	}
}

std::string convertIntToString(int a) {
	if (a == 0) {
		std::string result = "0";
		return result;
	}
	if (a > 0) {
		std::string result = "";
		while (a > 0) {
			int last = a % 10;
			a = a / 10;
			char c = last + '0';
			result = c + result;
		}
		return result;
	}
	a = -a;
	std::string result = "-";
	return result + convertIntToString(a);
}