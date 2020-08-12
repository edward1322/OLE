#ifndef UTILITIES_H
#define UTILITIES_H


#include <complex>

class utilities
{
public:
    utilities();
    ~utilities();


    double ComplexAbs(std::complex<double> a);
    double ComplexAbsSquared(std::complex<double> a);
    double SimpsonsWeight (int i, int n);
};

#endif // UTILITIES_H
