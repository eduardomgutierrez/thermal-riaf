#include "luminosityPhotoHadronic.h"

#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
//#include <fmath/interpolation.h>
#include <flosses/crossSectionInel.h>
#include <finjection/pgammaPionInj.h>
#include <fmath/physics.h> 

//double lossesPhotoHadronic(double E, Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& psc, double phEmin, double phEmax)
//{  //E=Ep

double luminosityPhotoHadronic(double Eg, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
								double phEmin, double phEmax)
{
	double diezEg = 10.0*Eg;
	double distCreator = 0.0;
	if (diezEg < p.emin() || diezEg> p.emax())
		distCreator = 0.0;
	else
		distCreator = p.distribution.interpolate({ { 0, diezEg } }, &psc); 

	double t_1 = t_pion_PHsimple(diezEg, p, tpf, psc, phEmin, phEmax);
	//double t_1 = t_pion_PHsimple(diezE,p,tpf,psc,phEmin,phEmax);
						//[&tpf, &psc](double x) {return tpf.interpolate({ {0, x } }, &psc); },
						//[&tpf, &phEmin,&phEmax, &psc](double x) {if (x < phEmin || x> phEmax){return 0.0;}
						//						else{ return tpf.interpolate({ {0, x } }, &psc);}},  
						//phEmin, phEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	double omega = omegaPHsimple(diezEg, p, tpf, psc, phEmin, phEmax);
	//double omega = omegaPHsimple(diezE,p,tpf,psc,phEmin,phEmax);
						//[&tpf, &psc](double x) {return tpf.interpolate({ {0, x } }, &psc); },
						//[&tpf, &phEmin,&phEmax, &psc](double x) {if (x < phEmin || x> phEmax){return 0.0;}
						//						else{ return tpf.interpolate({ {0, x } }, &psc);}}, 
						//phEmin, phEmax);
	if (omega > 0.0 && t_1 > 0.0) {
		double averageInel = t_1/omega;
		double k1 = 0.2;
		double k2 = 0.6;
		double p1 = (k2-averageInel)/(k2-k1);
        double p2 = 1.0 - p1;
        double nChargedPion = 0.5*p1 + 2.*p2;
		double nNeutralPion = (1.0-0.5)*p1 + p2;
        
        double luminosity = 20.0 * nNeutralPion * distCreator * omega;
        
        //double lumPion = P2(5.) * (t_1 * nNeutralPion/(nChargedPion+nNeutralPion)) * distCreator;
		//double luminosity = P2(2.) * lumPion;
		return luminosity*P2(Eg);
	} else {
		return 0.0;
    }
}


/*double luminosityPhotoHadronic2(double Eg, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
								double phEmin, double phEmax)
{
	
    double integral = integSimpsonLog(p.emin(), p.emax(),
                        [&psc, &p, phEmin, phEmax, &tpf, Eg]
                        (double Ep)
                        {
                            double Np = p.distribution.interpolate({{DIM_E, Ep}}, &psc);
                            double eth = 
                            return Np/Ep * integ2;
                        }, 100);
}*/