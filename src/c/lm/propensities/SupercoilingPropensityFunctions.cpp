/*
 * University of Illinois Open Source License
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to
 * do so, subject to the following conditions:
 *
 * - Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimers.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the names of the Roberts Group, Johns Hopkins University,
 * nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts, Max Klein
 */

#include <list>
#include <map>
#include <string>
#include <vector>
#include "math.h"

#include "lm/ClassFactory.h"
#include "lm/me/PropensityFunction.h"
#include "lm/propensities/SupercoilingPropensityFunctions.h"


namespace lm {
namespace propensities {

bool SupercoilingPropensityFunctions::registered=SupercoilingPropensityFunctions::registerClass();

bool SupercoilingPropensityFunctions::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::me::PropensityFunctionCollection", "lm::propensities::SupercoilingPropensityFunctions", &SupercoilingPropensityFunctions::allocateObject);
    return true;
}

void* SupercoilingPropensityFunctions::allocateObject()
{
    return new SupercoilingPropensityFunctions();
}

SupercoilingPropensityFunctions::SupercoilingPropensityFunctions()
{
}

SupercoilingPropensityFunctions::~SupercoilingPropensityFunctions()
{
}

class PoissonStepPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10000;

    PoissonStepPropensity(uint s1, uint s2, double k, int n) :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),k(k),n(n) {}
    uint s1,s2;
    double k;
    int n;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        int x1 = speciesCounts[s1];
        int x2 = speciesCounts[s2];
        return (x2 >= n)?(k*static_cast<double>(x1)):(0.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "poisson step propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "poisson step propensity require one species two dependency, had",s2.len);

        // Find the rate costant.
        if (k.len != 2)  throw InvalidArgException("k", "poisson step propensity needs two parameters, had",k.len);
        if (k[1] < 1.0)  throw InvalidArgException("k[1]", "poisson step propensity needs at least one step, had",k[1]);

        return new PoissonStepPropensity(s1[0], s2[0], k[0], static_cast<int>(k[1]));
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {
            "k1 * x1 * x2^k2",
        NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "PoissonStepPropensity", expressions, &create);
    }
};

class RNAPStallPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10001;

    RNAPStallPropensity(uint s1, uint s2, uint s3, double lk0, double torqueDomainBoundary0, double torqueDomainBoundary1, double torqueValue0, double torqueValue1, double torqueValue2, double stallCutoff, double stallRate, double ncoilCutoff, double pcoilCutoff)
    :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),s3(s3),lk0(lk0), torqueDomainBoundary0(torqueDomainBoundary0),torqueDomainBoundary1(torqueDomainBoundary1),torqueValue0(torqueValue0),torqueValue1(torqueValue1),torqueValue2(torqueValue2),stallCutoff(stallCutoff), stallRate(stallRate), ncoilCutoff(ncoilCutoff), pcoilCutoff(pcoilCutoff)
    {}

    uint s1, s2, s3;
    double lk0, torqueDomainBoundary0, torqueDomainBoundary1, torqueValue0, torqueValue1, torqueValue2, stallCutoff, stallRate, ncoilCutoff, pcoilCutoff;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        int RNAP = speciesCounts[s1];
		int Lk_up = speciesCounts[s2];
		int Lk_down = speciesCounts[s3];

		double sigma_up = (Lk_up - lk0)/lk0;
		double sigma_dn = (Lk_down - lk0)/lk0;
	
		double tau_up_t = piecewiseTorque(sigma_up);
		double tau_up = (tau_up_t<-10.5?-10.5:tau_up_t);
		double tau_dn = piecewiseTorque(sigma_dn);
		double tau = tau_dn - tau_up;

		if (sigma_up <= ncoilCutoff)
			return RNAP*stallRate;
		else if (sigma_dn >= pcoilCutoff)
			return RNAP*stallRate;
		else
			return RNAP*(tau>=stallCutoff?stallRate:0.0);
    }

    double piecewiseTorque(double sigma) const
    {
        if (fabs(sigma) < torqueDomainBoundary0)
            return torqueValue0*sigma;
        else if (fabs(sigma) <= torqueDomainBoundary1)
            return torqueValue1*(sigma<0?-1.0:1.0);
        else
            return torqueValue2*sigma;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "poisson step propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "poisson step propensity require one species two dependency, had",s2.len);
        utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
        if (s3.len != 1) throw InvalidArgException("D", "poisson step propensity require one species three dependency, had",s3.len);

        // Find the rate costant.
        if (k.len != 10)  throw InvalidArgException("k", "RNAP stall propensity needs seven parameters, had",k.len);

        return new RNAPStallPropensity(s1[0], s2[0], s3[0], k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "RNAPStallPropensity", expressions, &create);
    }
};


class RNAPResumePropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10002;

    RNAPResumePropensity(uint s1, uint s2, uint s3, double lk0, double torqueDomainBoundary0, double torqueDomainBoundary1, double torqueValue0, double torqueValue1, double torqueValue2, double stallCutoff, double resumeRate, double pcoilCutoff, double ncoilCutoff)
    :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),s3(s3),lk0(lk0),torqueDomainBoundary0(torqueDomainBoundary0),torqueDomainBoundary1(torqueDomainBoundary1),torqueValue0(torqueValue0),torqueValue1(torqueValue1),torqueValue2(torqueValue2),stallCutoff(stallCutoff),resumeRate(resumeRate),pcoilCutoff(pcoilCutoff),ncoilCutoff(ncoilCutoff)
    {}

    uint s1, s2, s3;
    double lk0, torqueDomainBoundary0, torqueDomainBoundary1, torqueValue0, torqueValue1, torqueValue2, stallCutoff, resumeRate, pcoilCutoff, ncoilCutoff;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        int RNAP = speciesCounts[s1];
		int Lk_up = speciesCounts[s2];
		int Lk_down = speciesCounts[s3];

		double sigma_up = (Lk_up - lk0)/lk0;
		double sigma_dn = (Lk_down - lk0)/lk0;

        double tau_up_t = piecewiseTorque(sigma_up);
		double tau_up = (tau_up_t<-10.5?-10.5:tau_up_t);
        double tau_dn = piecewiseTorque(sigma_dn);
        double tau = tau_dn - tau_up;

		if (sigma_up <= ncoilCutoff)
			return 0.0;
		else if (sigma_dn >= pcoilCutoff)
			return 0.0;
		else
			return RNAP*(tau>=stallCutoff?0.0:resumeRate);
//		double up_index = 1.0*(sigma_up>ncoilCutoff?1.0:0.0);
//		double down_index = 1.0*(sigma_dn<pcoilCutoff?1.0:0.0);
//       return up_index*down_index*RNAP*(tau<stallCutoff?resumeRate:0.0);
    }

    double piecewiseTorque(double sigma) const
    {
        if (fabs(sigma) < torqueDomainBoundary0)
            return torqueValue0*sigma;
        else if (fabs(sigma) <= torqueDomainBoundary1)
            return torqueValue1*(sigma<0?-1.0:1.0);
        else
            return torqueValue2*sigma;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "poisson step propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "poisson step propensity require one species two dependency, had",s2.len);
        utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
        if (s3.len != 1) throw InvalidArgException("D", "poisson step propensity require one species three dependency, had",s3.len);

        // Find the rate costant.
        if (k.len != 10)  throw InvalidArgException("k", "RNAP resume propensity needs seven parameters, had",k.len);

        return new RNAPResumePropensity(s1[0], s2[0], s3[0], k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "RNAPResumePropensity", expressions, &create);
    }
};



class InitiationPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10003;
    InitiationPropensity(uint s1, uint s2, uint s3, double Boundary0, double Boundary1, double initiationRate, double lk0, double basalLevel)
	:PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),s3(s3),Boundary0(Boundary0), Boundary1(Boundary1), initiationRate(initiationRate), lk0(lk0), basalLevel(basalLevel)
	{}
	
    uint s1, s2, s3;
	double Boundary0, Boundary1, initiationRate, lk0, basalLevel;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int DNA = speciesCounts[s1];
		int Lk = speciesCounts[s2];
		int DNA_near = speciesCounts[s3];

		double Sigma = (Lk-lk0)/lk0;

		if (Sigma>=Boundary0)
			return DNA*initiationRate*basalLevel*DNA_near;
		else if (Sigma<=Boundary1)
			return DNA*initiationRate*DNA_near;
		else
			return DNA*initiationRate*((Boundary0-Sigma)*(1-basalLevel)/(Boundary0-Boundary1) + basalLevel)*DNA_near;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "initiation propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "initiation propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "initiation propensity require one species two dependency, had",s3.len);
        // Find the rate costant.
        if (k.len != 5)  throw InvalidArgException("k", "initiation propensity needs three parameters, had",k.len);

        return new InitiationPropensity(s1[0], s2[0], s3[0], k[0], k[1], k[2], k[3], k[4]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "InitiationPropensity", expressions, &create);
    }
};



class GyraseCatalysisPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10004;

	GyraseCatalysisPropensity(uint s1, uint s2, double sigma_t, double epsilon, double catalyticRate, double lk0)
    :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),sigma_t(sigma_t),epsilon(epsilon),catalyticRate(catalyticRate), lk0(lk0)
    {}

    uint s1, s2;
	double sigma_t, epsilon, catalyticRate, lk0;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int Gyrase_bind = speciesCounts[s1];
		int Lk = speciesCounts[s2];

		double Sigma = (Lk-lk0)/lk0;

		return catalyticRate*Gyrase_bind/(1+exp (-(Sigma-sigma_t)/epsilon))*(Lk>0?1.0:0.0);
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "gyrase catalytic propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "gyrase catalytic propensity require one species two dependency, had",s2.len);

        // Find the rate costant.
        if (k.len != 4)  throw InvalidArgException("k", "gyrase catalytic propensity needs three parameters, had",k.len);

        return new GyraseCatalysisPropensity(s1[0], s2[0], k[0], k[1], k[2], k[3]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "GyraseCatalysisPropensity", expressions, &create);
    }
};




class SCDiffusionPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10005;

	SCDiffusionPropensity(uint s1, uint s2, double diffusionCoefficient)
	:PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2), diffusionCoefficient(diffusionCoefficient)
    {}

    uint s1, s2;
	double diffusionCoefficient;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int SC = speciesCounts[s1];
		int DNA = speciesCounts[s2];
		return DNA*(SC>0?1.0:0.0)*diffusionCoefficient;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "SC diffusion propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "SC diffusion propensity require one species two dependency, had",s2.len);

        // Find the rate costant.
        if (k.len != 1)  throw InvalidArgException("k", "SC diffusion propensity needs three parameters, had",k.len);

        return new SCDiffusionPropensity(s1[0], s2[0], k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "SCDiffusionPropensity", expressions, &create);
    }
};



class TopoICatalysisPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10006;

	TopoICatalysisPropensity(uint s1, uint s2, double catalyticRate, double lk0)
    :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),catalyticRate(catalyticRate), lk0(lk0)
    {}

    uint s1, s2;
	double catalyticRate, lk0;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int TopoI_bind = speciesCounts[s1];
		int Lk = speciesCounts[s2];

		return TopoI_bind*(Lk>lk0?0.0:1.0)*catalyticRate;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "TopoI catalytic propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "TopoI catalytic propensity require one species two dependency, had",s2.len);

        // Find the rate costant.
        if (k.len != 2)  throw InvalidArgException("k", "TopoI catalytic propensity needs three parameters, had",k.len);

        return new TopoICatalysisPropensity(s1[0], s2[0], k[0], k[1]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TopoICatalysisPropensity", expressions, &create);
    }
};


class SCDiffusionPropensityNew : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10007;

	SCDiffusionPropensityNew(uint s1, uint s2, uint s3, double diffusionCoefficient)
	:PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2), s3(s3),diffusionCoefficient(diffusionCoefficient)
    {}

    uint s1, s2, s3;
	double diffusionCoefficient;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int SC = speciesCounts[s1];
		int SC_near = speciesCounts[s2];
		int DNA = speciesCounts[s3];

		return DNA*(SC-SC_near)*(SC>SC_near?1.0:0.0)*diffusionCoefficient;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "SC diffusion propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "SC diffusion propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "SC diffusion propensity require one species three dependency, had",s3.len);

        // Find the rate costant.
        if (k.len != 1)  throw InvalidArgException("k", "SC diffusion propensity needs three parameters, had",k.len);

        return new SCDiffusionPropensityNew(s1[0], s2[0], s3[0], k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "SCDiffusionPropensityNew", expressions, &create);
    }
};


class RNAPtranslocation : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10008;

	RNAPtranslocation(uint s1, uint s2, uint s3, double translocationRate)
	:PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2), s3(s3), translocationRate(translocationRate)
    {}

    uint s1, s2, s3;
	double translocationRate;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int DNA_near = speciesCounts[s1];
		int RNAP_tmp = speciesCounts[s2];
		int DNA_near2 = speciesCounts[s3];

		return DNA_near*RNAP_tmp*(DNA_near2>0?1.0:0.0)*translocationRate;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "RNAP translocation propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "RNAP translocation propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "RNAP translocation propensity require one species three dependency, had",s3.len);

        // Find the rate costant.
        if (k.len != 1)  throw InvalidArgException("k", "RNAP translocation propensity needs three parameters, had",k.len);

        return new RNAPtranslocation(s1[0], s2[0], s3[0], k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "RNAPtranslocation", expressions, &create);
    }
};



class TopoIbinding : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10009;

	TopoIbinding(uint s1, uint s2, uint s3, uint s4, uint s5, uint s6, uint s7, uint s8, double topoI_binding_rate1, double topoI_binding_rate2)
	:PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2), s3(s3), s4(s4), s5(s5),s6(s6), s7(s7), s8(s8), topoI_binding_rate1(topoI_binding_rate1), topoI_binding_rate2(topoI_binding_rate2)
    {}

    uint s1, s2, s3, s4, s5, s6, s7, s8;
	double topoI_binding_rate1;
	double topoI_binding_rate2;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int TopoI_unbind = speciesCounts[s1];
		int DNA = speciesCounts[s2];
		int RNAP_up = speciesCounts[s3];
		int RNAP_tmp_up = speciesCounts[s4];
		int RNAP_stall_up = speciesCounts[s5];
		int RNAP_down = speciesCounts[s6];
		int RNAP_tmp_down = speciesCounts[s7];
		int RNAP_stall_down = speciesCounts[s8];	
	
		if (RNAP_up+RNAP_tmp_up+RNAP_stall_up+RNAP_down+RNAP_tmp_down+RNAP_stall_down==0)
			return TopoI_unbind*DNA*topoI_binding_rate1;
		else
			return TopoI_unbind*DNA*topoI_binding_rate2;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species three dependency, had",s3.len);
		utuple s4 = getSpecificDependencies(reactionIndex, D, 4);
		if (s4.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species three dependency, had",s4.len);
		utuple s5 = getSpecificDependencies(reactionIndex, D, 5);
		if (s5.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species three dependency, had",s5.len);
		utuple s6 = getSpecificDependencies(reactionIndex, D, 6);
		if (s6.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species three dependency, had",s6.len);
		utuple s7 = getSpecificDependencies(reactionIndex, D, 7);
		if (s7.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species three dependency, had",s7.len);
		utuple s8 = getSpecificDependencies(reactionIndex, D, 8);
		if (s8.len != 1) throw InvalidArgException("D", "TopoI binding propensity require one species three dependency, had",s8.len);
        // Find the rate costant.
        if (k.len != 2)  throw InvalidArgException("k", "TopoI binding propensity needs three parameters, had",k.len);

        return new TopoIbinding(s1[0], s2[0], s3[0], s4[0], s5[0], s6[0], s7[0], s8[0], k[0], k[1]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TopoIbinding", expressions, &create);
    }
};



class mRNA_degradation_elongation : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10010;

	mRNA_degradation_elongation(uint s1, uint s2, uint s3, double RNase_elongation_rate)
	:PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2), s3(s3), RNase_elongation_rate(RNase_elongation_rate)
    {}

    uint s1, s2, s3;
	double RNase_elongation_rate;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int mRNA_k = speciesCounts[s1];
		int mRNA_degr_k = speciesCounts[s2];
		int mRNA_pre_k = speciesCounts[s3];
		
		if (mRNA_pre_k==0)
			return mRNA_k*mRNA_degr_k*RNase_elongation_rate;
		else
			return 0;
   }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "mRNA_degradation_elongation propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "mRNA_degradation_elongation propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "mRNA_degradation_elongation propensity require one species three dependency, had",s3.len);
        // Find the rate costant.
        if (k.len != 1)  throw InvalidArgException("k", "mRNA_degradation_elongation propensity needs three parameters, had",k.len);

        return new mRNA_degradation_elongation(s1[0], s2[0], s3[0], k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "mRNA_degradation_elongation", expressions, &create);
    }
};


class ClearPropensity : public lm::me::PropensityFunction

{
public:
    static const uint REACTION_TYPE = 10011;

	ClearPropensity(uint s1, uint s2, double clearRate)
    :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),clearRate(clearRate)
    {}

    uint s1, s2;
	double clearRate;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int clear = speciesCounts[s1];
		int stall_start = speciesCounts[s2];

		if (stall_start>0)
			return 0;
		else
			return clear*clearRate;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "clear propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "clear propensity require one species two dependency, had",s2.len);
        // Find the rate costant.
        if (k.len != 1)  throw InvalidArgException("k", "clear propensity needs one parameter, had",k.len);

        return new ClearPropensity(s1[0], s2[0], k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "ClearPropensity", expressions, &create);
    }
};


class PrematurePropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10012;

	PrematurePropensity(uint s1, uint s2, uint s3, double prematureThresh, double prematureRate, double lk0, double torqueDomainBoundary0, double torqueDomainBoundary1, double torqueValue0, double torqueValue1, double torqueValue2)
    :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),s3(s3),prematureThresh(prematureThresh),prematureRate(prematureRate),lk0(lk0), torqueDomainBoundary0(torqueDomainBoundary0),torqueDomainBoundary1(torqueDomainBoundary1),torqueValue0(torqueValue0),torqueValue1(torqueValue1),torqueValue2(torqueValue2)
    {}

    uint s1, s2, s3;
	double prematureThresh, prematureRate, lk0, torqueDomainBoundary0, torqueDomainBoundary1, torqueValue0, torqueValue1, torqueValue2;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int RNAP_stall = speciesCounts[s1];
		int Lk_up = speciesCounts[s2];
		int Lk_down = speciesCounts[s3];

		double sigma_up = (Lk_up - lk0)/lk0;
		
		return RNAP_stall*(sigma_up<prematureThresh?prematureRate:0.0);
	}

	double piecewiseTorque(double sigma) const
    {
        if (fabs(sigma) < torqueDomainBoundary0)
            return torqueValue0*sigma;
        else if (fabs(sigma) <= torqueDomainBoundary1)
            return torqueValue1*(sigma<0?-1.0:1.0);
        else
            return torqueValue2*sigma;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "clear propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "clear propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "clear propensity require one species two dependency, had",s3.len);
        // Find the rate costant.
        if (k.len != 8)  throw InvalidArgException("k", "clear propensity needs one parameter, had",k.len);

        return new PrematurePropensity(s1[0], s2[0], s3[0], k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "PrematurePropensity", expressions, &create);
    }
};



class PrematureResumePropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 10013;

	PrematureResumePropensity(uint s1, uint s2, uint s3, double prematureThresh, double prematureRate, double lk0, double torqueDomainBoundary0, double torqueDomainBoundary1, double torqueValue0, double torqueValue1, double torqueValue2)
    :PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),s3(s3),prematureThresh(prematureThresh),prematureRate(prematureRate),lk0(lk0), torqueDomainBoundary0(torqueDomainBoundary0),torqueDomainBoundary1(torqueDomainBoundary1),torqueValue0(torqueValue0),torqueValue1(torqueValue1),torqueValue2(torqueValue2)
    {}

    uint s1, s2, s3;
	double prematureThresh, prematureRate, lk0, torqueDomainBoundary0, torqueDomainBoundary1, torqueValue0, torqueValue1, torqueValue2;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int RNAP_stall = speciesCounts[s1];
		int Lk_up = speciesCounts[s2];
		int Lk_down = speciesCounts[s3];

		double sigma_up = (Lk_up - lk0)/lk0;
		
		return RNAP_stall*(sigma_up>prematureThresh?prematureRate:0.0);
	}

	double piecewiseTorque(double sigma) const
    {
        if (fabs(sigma) < torqueDomainBoundary0)
            return torqueValue0*sigma;
        else if (fabs(sigma) <= torqueDomainBoundary1)
            return torqueValue1*(sigma<0?-1.0:1.0);
        else
            return torqueValue2*sigma;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "clear propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "clear propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "clear propensity require one species two dependency, had",s3.len);
        // Find the rate costant.
        if (k.len != 8)  throw InvalidArgException("k", "clear propensity needs one parameter, had",k.len);

        return new PrematureResumePropensity(s1[0], s2[0], s3[0], k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "PrematureResumePropensity", expressions, &create);
    }
};



class PIPropensity : public lm::me::PropensityFunction

{
public:
    static const uint REACTION_TYPE = 10014;

	PIPropensity(uint s1, uint s2, uint s3, double inhibitionRate)
	:PropensityFunction(REACTION_TYPE,0),s1(s1),s2(s2),s3(s3),inhibitionRate(inhibitionRate)
    {}

    uint s1, s2, s3;
	double inhibitionRate;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
		int RNAP_index = speciesCounts[s1];
		int DNA_promoter = speciesCounts[s2];
		int DNA_promoter_dn = speciesCounts[s3];

		return RNAP_index*DNA_promoter*DNA_promoter_dn*inhibitionRate;
    }


#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple s1 = getSpecificDependencies(reactionIndex, D, 1);
        if (s1.len != 1) throw InvalidArgException("D", "PI propensity require one species one dependency, had",s1.len);
        utuple s2 = getSpecificDependencies(reactionIndex, D, 2);
        if (s2.len != 1) throw InvalidArgException("D", "PI propensity require one species two dependency, had",s2.len);
		utuple s3 = getSpecificDependencies(reactionIndex, D, 3);
		if (s3.len != 1) throw InvalidArgException("D", "PI propensity require one species two dependency, had",s3.len);
        // Find the rate costant.
        if (k.len != 1)  throw InvalidArgException("k", "PI propensity needs one parameter, had",k.len);

        return new PIPropensity(s1[0], s2[0], s3[0], k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "PIPropensity", expressions, &create);
    }
};


list<lm::me::PropensityFunctionDefinition> SupercoilingPropensityFunctions::getPropensityFunctionDefinitions()
{
    list<lm::me::PropensityFunctionDefinition> defs;
    defs.push_back(PoissonStepPropensity::registerFunction());
    defs.push_back(RNAPStallPropensity::registerFunction());
    defs.push_back(RNAPResumePropensity::registerFunction());
	defs.push_back(InitiationPropensity::registerFunction());
	defs.push_back(GyraseCatalysisPropensity::registerFunction());
	defs.push_back(SCDiffusionPropensity::registerFunction());
	defs.push_back(TopoICatalysisPropensity::registerFunction());
	defs.push_back(SCDiffusionPropensityNew::registerFunction());
	defs.push_back(RNAPtranslocation::registerFunction());
	defs.push_back(TopoIbinding::registerFunction());
	defs.push_back(mRNA_degradation_elongation::registerFunction());
	defs.push_back(ClearPropensity::registerFunction());
//	defs.push_back(PrematurePropensity::registerFunction());
//	defs.push_back(PrematureResumePropensity::registerFunction());
	defs.push_back(PIPropensity::registerFunction());
	return defs;
}

}
}
