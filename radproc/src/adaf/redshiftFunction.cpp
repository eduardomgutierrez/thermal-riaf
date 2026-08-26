#include <iostream>
#include <cmath>

#include "redshiftFunction.h"
#include "adafFunctions.h"

extern double schwRadius;
extern double gravRadius;
extern double blackHoleSpin;
extern double horizonRadius;
extern int readPrecomputedADAF;
extern size_t nR, nRcd;

Vector redshift_to_inf;							// Vector [jR] with the redshift between a shell and infinity.
Vector redshift_CD_to_inf;						// Vector [jRcd] with the redshift between the CD and infinity.
Matrix redshift;								// Matrix [jR][j'R] with the redshift between shells.
Matrix redshift_CD_to_RIAF;						// Matrix [jRcd][jR] with the redshift between the CD and a shell
												// of the RIAF.
Matrix redshift_RIAF_to_CD;						// Matrix [jR][jRcd] with the redshift between a shell of the
												// RIAF and the CD.

// Compute redshift factor for Kerr metric
// redshift_factor = Delta^(1/2) * r / (gamma_r * gamma_phi * A^(1/2))
// where Delta = r^2 - 2r + a^2, A = (r^2 + a^2)^2 - a^2*Delta
// gamma_phi = 1/sqrt(1 - (v_phi)^2), gamma_r = 1/sqrt(1 - V^2)
// V = v_r / gamma_phi, v_r is radial velocity (LNRF), v_phi is azimuthal velocity (ZAMO)
double kerrRedshiftFactor(double r_cm)
{
	double a = blackHoleSpin;
	double r = r_cm / gravRadius;  // Convert to dimensionless r/r_g
	
	// Kerr metric functions
	double r2 = r * r;
	double a2 = a * a;
	double Delta = r2 - 2.0*r + a2;
	double A = (r2 + a2)*(r2 + a2) - a2*Delta;
	
	if (Delta <= 0.0 || A <= 0.0) return 1.0;  // Inside or at horizon
	
	// Get radial velocity v^r (in units of c, as seen by LNRF)
	double vr = radialVel(r_cm) / cLight;  // radialVel returns v in cm/s
	if (std::abs(vr) >= 1.0) vr = -0.99;
	
	// Use azimuthal velocity profile directly from precomputed ADAF data.
	double vphi = precompAzimuthalVel(r_cm) / cLight;  // in units of c
	
	if (std::abs(vphi) >= 1.0) vphi = 0.99 * (vphi > 0 ? 1.0 : -1.0);
	
	// Lorentz factors
	double gamma_phi2 = 1.0 / (1.0 - vphi*vphi);
	double gamma_phi = sqrt(gamma_phi2);
	double V = vr; // / gamma_phi;  // Effective radial velocity
	if (std::abs(V) >= 1.0) V = 0.99 * (V < 0 ? -1.0 : 1.0);  // Safety check
	double gamma_r = 1.0 / sqrt(1.0 - V*V);
	
	// Redshift factor
	double redshift_factor = sqrt(Delta) * r / (gamma_r * gamma_phi * sqrt(A));

  std::cout << "Kerr vphi at r = " << r << " (in units of r_g): " << vphi << std::endl;
  std::cout << "Kerr redshift factor at r = " << r << " (in units of r_g): " << redshift_factor << std::endl;
	
	return (redshift_factor > 0.0 && std::isfinite(redshift_factor)) ? redshift_factor : 1.0;
}

// Compute Schwarzschild/Paczynskii-Wiita redshift factor (original formula)
double schwarzschildRedshiftFactor(double r_cm)
{
	double vr = radialVel(r_cm);
	double beta = vr / cLight;
	if (std::abs(beta) >= 1.0) beta = -0.9;
	double redshift_factor = sqrt((1.0 - schwRadius/r_cm) * (1.0 - beta*beta));
	return (redshift_factor > 0.0) ? redshift_factor : 1.0;
}

void redshiftFactor(State& st)
{
	redshift_to_inf.resize(nR,1.0);
	redshift_CD_to_inf.resize(nRcd,1.0);
	matrixInit(redshift,nR,nR,1.0);
	matrixInit(redshift_CD_to_RIAF,nRcd,nR,1.0);
	matrixInit(redshift_RIAF_to_CD,nR,nRcd,1.0);
	
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R);
		double vr = radialVel(r);
		double beta = vr/cLight;
		if (abs(beta) >= 1.0) beta = -0.99;
		
		// Use Kerr formula for precomputed ADAF, otherwise use Schwarzschild/Paczynskii-Wiita
		double redshift_factor;
		if (readPrecomputedADAF && blackHoleSpin != 0.0) {
			redshift_factor = kerrRedshiftFactor(r);
		} else {
			redshift_factor = schwarzschildRedshiftFactor(r);
		}
		redshift_to_inf[jR] = redshift_factor;
		cout << "r [r_g] = " << r/gravRadius << "\t (1+z)^-1 = " << redshift_to_inf[jR] << endl;
		size_t jjR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itRR) {
			double rr = itRR.val(DIM_R);
			double vrr = radialVel(rr);
			double denom = 1.0 - vr*vrr/cLight2;
			if (std::abs(denom) < 1e-10) denom = 1e-10;  // Safety check
			double relative_velocity = abs(vr-vrr) / denom;
			double relative_beta = relative_velocity/cLight;
			if (abs(relative_beta) >= 1.0) relative_beta = 0.9;
			double doppler_factor = sqrt( (1.0-relative_beta) / (1.0+relative_beta) );
			
			// Gravitational redshift factor between shells
			double grav_factor_shells;
			if (readPrecomputedADAF && blackHoleSpin != 0.0) {
				// For Kerr: use ratio of redshift factors (gravitational part)
				double a = blackHoleSpin;
				double r_rg = r / gravRadius;
				double rr_rg = rr / gravRadius;
				double Delta_r = r_rg*r_rg - 2.0*r_rg + a*a;
				double Delta_rr = rr_rg*rr_rg - 2.0*rr_rg + a*a;
				double A_r = pow(r_rg*r_rg + a*a, 2) - a*a*Delta_r;
				double A_rr = pow(rr_rg*rr_rg + a*a, 2) - a*a*Delta_rr;
				if (Delta_r > 0 && Delta_rr > 0 && A_r > 0 && A_rr > 0) {
					grav_factor_shells = sqrt(Delta_r/Delta_rr) * (r_rg/rr_rg) * sqrt(A_rr/A_r);
				} else {
					grav_factor_shells = 1.0;
				}
			} else {
				double arg1 = 1.0-schwRadius/r;
				double arg2 = 1.0-schwRadius/rr;
				if (arg1 <= 0.0) arg1 = 1e-6;
				if (arg2 <= 0.0) arg2 = 1e-6;
				grav_factor_shells = sqrt( arg1 / arg2 );
			}
			double redshift_val = grav_factor_shells * doppler_factor;
			// Safety check: redshift must be positive and finite
			if (!std::isfinite(redshift_val) || redshift_val <= 0.0) {
				redshift_val = 1.0;
			}
			redshift[jR][jjR] = redshift_val;
			jjR++;
		},{0,-1,0});
		
		size_t jRcd=0;
		st.photon.ps.iterate([&](const SpaceIterator& itRcd) {
			double rCD = itRcd.val(DIM_Rcd);
			double beta_local = beta;
			if (rCD < r) beta_local = -beta_local;
			double doppler_factor = sqrt( (1.0-beta_local) / (1.0+beta_local) );
			
			// Gravitational redshift factor between RIAF shell and CD
			double redshift_factor_CD;
			if (readPrecomputedADAF && blackHoleSpin != 0.0) {
				double a = blackHoleSpin;
				double r_rg = r / gravRadius;
				double rCD_rg = rCD / gravRadius;
				double Delta_r = r_rg*r_rg - 2.0*r_rg + a*a;
				double Delta_rCD = rCD_rg*rCD_rg - 2.0*rCD_rg + a*a;
				double A_r = pow(r_rg*r_rg + a*a, 2) - a*a*Delta_r;
				double A_rCD = pow(rCD_rg*rCD_rg + a*a, 2) - a*a*Delta_rCD;
				if (Delta_r > 0 && Delta_rCD > 0 && A_r > 0 && A_rCD > 0) {
					redshift_factor_CD = sqrt(Delta_r/Delta_rCD) * (r_rg/rCD_rg) * sqrt(A_rCD/A_r);
				} else {
					redshift_factor_CD = 1.0;
				}
			} else {
				double arg1 = 1.0-schwRadius/r;
				double arg2 = 1.0-schwRadius/rCD;
				if (arg1 <= 0.0) arg1 = 1e-6;
				if (arg2 <= 0.0) arg2 = 1e-6;
				redshift_factor_CD = sqrt( arg1 / arg2 );
			}
			double redshift_val_CD = redshift_factor_CD * doppler_factor;
			// Safety check: redshift must be positive and finite
			if (!std::isfinite(redshift_val_CD) || redshift_val_CD <= 0.0) {
				redshift_val_CD = 1.0;
			}
			redshift_RIAF_to_CD[jR][jRcd] = redshift_val_CD;
			jRcd++;
		},{0,0,-1});
		
		jR++;
	},{0,-1,0});
	
	size_t jRcd=0;
	st.photon.ps.iterate([&](const SpaceIterator& itRcd) {
		double rCD = itRcd.val(DIM_Rcd);
		double betaCD = rCD*keplAngVel(rCD)/cLight;
		if (std::abs(betaCD) >= 1.0) betaCD = 0.9 * (betaCD > 0 ? 1.0 : -1.0);  // Safety check
		
		// Redshift factor for CD to infinity
		double redshift_factor_CD_to_inf;
		if (readPrecomputedADAF && blackHoleSpin != 0.0) {
			double a = blackHoleSpin;
			double rCD_rg = rCD / gravRadius;
			double Delta_rCD = rCD_rg*rCD_rg - 2.0*rCD_rg + a*a;
			double A_rCD = pow(rCD_rg*rCD_rg + a*a, 2) - a*a*Delta_rCD;
			if (Delta_rCD > 0 && A_rCD > 0) {
				// For Keplerian disk, approximate gamma_phi from betaCD
				double gamma_phi = 1.0 / sqrt(1.0 - betaCD*betaCD);
				redshift_factor_CD_to_inf = sqrt(Delta_rCD) * rCD_rg / (gamma_phi * sqrt(A_rCD));
			} else {
				redshift_factor_CD_to_inf = 1.0;
			}
		} else {
			double arg1 = 1.0-schwRadius/rCD;
			double arg2 = 1.0-betaCD*betaCD;
			if (arg1 <= 0.0) arg1 = 1e-6;
			if (arg2 <= 0.0) arg2 = 1e-6;
			redshift_factor_CD_to_inf = sqrt( arg1 * arg2 );
		}
		// Safety check: redshift must be positive and finite
		if (!std::isfinite(redshift_factor_CD_to_inf) || redshift_factor_CD_to_inf <= 0.0) {
			redshift_factor_CD_to_inf = 1.0;
		}
		redshift_CD_to_inf[jRcd] = redshift_factor_CD_to_inf;
		
		size_t jjR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itRR) {
			double rr = itRR.val(DIM_R);
			double vrr = radialVel(rr);
			double beta = -vrr/cLight;
			if (abs(beta) >= 1.0) beta = 0.9;
			if (rCD < rr) beta = -beta;
			double doppler_factor = sqrt( (1.0-beta) / (1.0+beta) );
			
			// Gravitational redshift factor CD to RIAF
			double grav_factor;
			if (readPrecomputedADAF && blackHoleSpin != 0.0) {
				double a = blackHoleSpin;
				double rCD_rg = rCD / gravRadius;
				double rr_rg = rr / gravRadius;
				double Delta_rCD = rCD_rg*rCD_rg - 2.0*rCD_rg + a*a;
				double Delta_rr = rr_rg*rr_rg - 2.0*rr_rg + a*a;
				double A_rCD = pow(rCD_rg*rCD_rg + a*a, 2) - a*a*Delta_rCD;
				double A_rr = pow(rr_rg*rr_rg + a*a, 2) - a*a*Delta_rr;
				if (Delta_rCD > 0 && Delta_rr > 0 && A_rCD > 0 && A_rr > 0) {
					grav_factor = sqrt(Delta_rCD/Delta_rr) * (rCD_rg/rr_rg) * sqrt(A_rr/A_rCD);
				} else {
					grav_factor = 1.0;
				}
			} else {
				double arg1 = 1.0-schwRadius/rCD;
				double arg2 = 1.0-schwRadius/rr;
				if (arg1 <= 0.0) arg1 = 1e-6;
				if (arg2 <= 0.0) arg2 = 1e-6;
				grav_factor = sqrt( arg1 / arg2 );
			}
			double redshift_val_CD_RIAF = grav_factor * doppler_factor;
			// Safety check: redshift must be positive and finite
			if (!std::isfinite(redshift_val_CD_RIAF) || redshift_val_CD_RIAF <= 0.0) {
				redshift_val_CD_RIAF = 1.0;
			}
			redshift_CD_to_RIAF[jRcd][jjR] = redshift_val_CD_RIAF;
			jjR++;
		},{0,-1,itRcd.coord[DIM_Rcd]});
		
		jRcd++;
	},{0,0,-1});
	
}