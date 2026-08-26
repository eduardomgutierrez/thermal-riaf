#include "ppPionInj.h"

#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <flosses/crossSectionInel.h>	

double Fpi_SYBYLL(double Epi, double Ep)         //funcion a integrar   x=Eproton; E=Epion
{ 
	double L       = log(Ep/1.602); //el 1.6 son TeV en erg
	double ap      = 3.67 + L*(0.83 + L*0.075);
	double Bpi     = ap + 0.25;
	double r       = 2.6/sqrt(ap);
	double alpha   = 0.98/sqrt(ap);
	double x  	   = Epi/Ep;
	double xalpha = pow(x, alpha);
	double factor  = 1.0 - xalpha;

	return (factor > 0. ? 
				4.0*alpha*Bpi * xalpha/x * pow(factor/(1.0+r*xalpha*factor), 4) * 
				(1.0/factor + r*(1.0-2.0*xalpha)/(1.0+r*xalpha*factor))*sqrt(1.0-chargedPionMass*cLight2/(x*Ep))
				: 0.0);
}

double Fpi_QGSJET(double Epi, double Ep)
{
    double L = log(Ep/1.602);
    double Bpi = 5.58 + L * (0.78 + 0.1*L);
    double r = 3.1/pow(Bpi, 1.5);
    double alpha = 0.89 / (sqrt(Bpi)*(1.0-exp(-0.33*Bpi)));
    double x = Epi / Ep;
    double xalpha = pow(x, alpha);
    double factor = 1.0 - xalpha;
    return (factor > 0. ?
                4.0*alpha*Bpi*xalpha/x * pow(factor/pow(1.0+r*xalpha, 3), 4) * 
                (1.0/factor + 3.0*r/(1.0+r*xalpha)) * sqrt(1.0 - chargedPionMass*cLight2/(x*Ep)) : 0.0);
}

double ppPionInj(double Epi, const Particle& creator, const double density, const SpaceCoord& psc)
{
    double inf = std::max(Epi, creator.emin());
	double sup = creator.emax();
	double integralP = integSimpsonLog(inf, sup,
                        [Epi, &creator, &psc]
                        (double Ep)
						{
                            double Np = creator.distribution.interpolate({{0, Ep}}, &psc);
                            double sigmapp = crossSectionHadronic(Ep);
							return Np/Ep * sigmapp * Fpi_QGSJET(Epi, Ep);
						}, 100);
	return 2.*integralP*cLight*density; 
}