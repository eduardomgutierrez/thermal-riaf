#pragma once

#include <fmath/physics.h>
#include <fparticle/Particle.h>

/* Pion injection due to pp interaction (Kelner et al. 2006) */

double ppPionInj(double Epi, const Particle& creator,
	const double density, const SpaceCoord& psc);
	