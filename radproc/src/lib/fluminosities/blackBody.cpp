#include "blackBody.h"
#include <fmath/physics.h>
#include <cmath>

double bb_RJ(double frequency, double temp) {
	if (temp <= 0.0 || frequency <= 0.0) return 0.0;
	return 2.0*frequency*frequency*boltzmann*temp / cLight2;
}

double bb(double frequency, double temp) {
	// Safety checks
	if (temp <= 0.0 || frequency <= 0.0) return 0.0;
	
	double x = planck*frequency/(boltzmann*temp);
	// Avoid overflow in exp for large x
	if (x > 700.0) return 0.0;  // Essentially zero at this point
	
	double expVal = exp(x) - 1.0;
	if (expVal <= 0.0) return 0.0;  // Safety
	
	double result = 2.0*planck*frequency*frequency*frequency / cLight2 / expVal;
	return std::isfinite(result) ? result : 0.0;
}