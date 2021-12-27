/*
 * University of Illinois Open Source License
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 *               Johns Hopkins University
 * 		 http://biophysics.jhu.edu/roberts/
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
 * Author(s): Elijah Roberts
 */

#include <list>
#include <map>
#include <string>
#include <vector>
#include "math.h"

#include "lm/ClassFactory.h"
#include "lm/me/PropensityFunction.h"
#include "lm/propensities/EnzymePropensityFunctions.h"


namespace lm {
namespace propensities {

bool EnzymePropensityFunctions::registered=EnzymePropensityFunctions::registerClass();

bool EnzymePropensityFunctions::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::me::PropensityFunctionCollection", "lm::propensities::EnzymePropensityFunctions", &EnzymePropensityFunctions::allocateObject);
    return true;
}

void* EnzymePropensityFunctions::allocateObject()
{
    return new EnzymePropensityFunctions();
}

EnzymePropensityFunctions::EnzymePropensityFunctions()
{
}

EnzymePropensityFunctions::~EnzymePropensityFunctions()
{
}



class FirstOrderMichaelisMenten : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 4002;

    FirstOrderMichaelisMenten(uint s1, double Rm, double Km) :PropensityFunction(REACTION_TYPE,1),s1(s1),Rm(Rm),Km(Km) {}
    uint s1;
    double Rm;
    double Km;

    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        return (Rm*double(speciesCounts[s1]))/(Km + double(speciesCounts[s1]));
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
        utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
        if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMenten needs one species dependencies, had", sd1.len);

        // Find the rate constant.
        if (k.len != 2) throw InvalidArgException("k", "FirstOrderMichaelisMenten needs two parameters, had", k.len);

        return new FirstOrderMichaelisMenten(sd1[0], k[0], k[1]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
	const char* unitsForConstants[] = {"item/second", "item", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderMichaelisMenten", "(k1 * x1) / (k2 + x1)", unitsForConstants, &create);
    }
};


class FirstOrderMichaelisMentenU1 : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4003;

	FirstOrderMichaelisMentenU1(uint s1, uint s2, double k, double Rm, double Km) :PropensityFunction(REACTION_TYPE,1),s1(s1),s2(s2),k(k),Rm(Rm),Km(Km) {}
	uint s1;
	uint s2;
	double k;
	double Rm;
	double Km;

	void changeVolume(double volumeMultiplier) {}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		double u1 = 1/(1+(k*double(speciesCounts[s1])));
		double MM = (Rm*double(speciesCounts[s2]))/(Km + double(speciesCounts[s2]));
		return u1*MM;
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU1 needs one species dependencies, had", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU1 needs one species dependencies, had", sd2.len);

		// Find the rate constant. 
		if (k.len < 3) throw InvalidArgException("k", "FirstOrderMichaelisMentenU1 needs three parameters, had", k.len);

		return new FirstOrderMichaelisMentenU1(sd1[0], sd2[0], k[0], k[1], k[2]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"1/item", "item/second", "item", NULL};
		const char* expressions[] = {"(1 / (1 + k1 * x1)) * ((k2 * x2) / (k3 +  x2))", "k2 * x2 * (1 / (1 + k1 * x1)) / (k3 + x2)", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderMichaelisMentenU1", expressions, unitsForConstants, &create);
	}
};

class FirstOrderDoubleMichaelisMenten : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4004;

	FirstOrderDoubleMichaelisMenten(uint s1, uint s2, double k1, double k2, double k3) :PropensityFunction(REACTION_TYPE,1),s1(s1),s2(s2),k1(k1),k2(k2),k3(k3)	{}
	uint s1;
	uint s2;
	double k1;
	double k2;
	double k3;

	void changeVolume(double volumeMultiplier) {}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		return k1*(double(speciesCounts[s1])/(k2 + double(speciesCounts[s1])))*(double(speciesCounts[s2])/(k3 + double(speciesCounts[s2])));
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderDoubleMichaelisMenten needs one species dependencies, had", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "FirstOrderDoubleMichaelisMenten needs one species dependencies, had", sd2.len);

		// Find the rate constant.
		if (k.len < 3) throw InvalidArgException("k", "FirstOrderDoubleMichaelisMenten needs three parameters, had", k.len);

		return new FirstOrderDoubleMichaelisMenten(sd1[0], sd2[0], k[0], k[1], k[2]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"item/second", "item", "item", NULL};
		const char* expressions[] = {"((k1 * x1) / (k2 + x1)) * (x2 / (k3 + x2))", "k1 * (x2 / (k3 + x2)) * (x1 / (k2 + x1))", "(k1 * (x2 / (k3 + x2))) * (x1) / (k2 + x1)", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderDoubleMichaelisMenten", expressions, unitsForConstants, &create);
	}
};

class FirstOrderProductSubstrateDependent2Species : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4005;

	FirstOrderProductSubstrateDependent2Species(uint s1, uint s2, double Km, double Rm) :PropensityFunction(REACTION_TYPE,1),s1(s1),s2(s2),Km(Km),Rm(Rm) {}
	uint s1;
	uint s2;
	double Km;
	double Rm;

	void changeVolume(double volumeMultiplier) {}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		double u4 = double(speciesCounts[s1])/(Km + double(speciesCounts[s1]));
		return u4*Rm*(double(speciesCounts[s2]) - double(speciesCounts[s1]));
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderProductSubstrateDependent2Species needs one species dependencies, had", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "FirstOrderProductSubstrateDependent2Species needs one species dependencies, had", sd2.len);

		// Find the rate constant.
		if (k.len < 2) throw InvalidArgException("k", "FirstOrderProductSubstrateDependent2Species needs two parameters, had", k.len);

		return new FirstOrderProductSubstrateDependent2Species(sd1[0], sd2[0], k[0], k[1]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"item", "1/second", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderProductSubstrateDependent2Species", "(x1 / (k1 + x1)) * k2 * (x2 - x1)", unitsForConstants, &create);
	}
};

class FirstOrderProductSubstrateDependent3Species : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4006;

	FirstOrderProductSubstrateDependent3Species(uint s1, uint s2, uint s3, double Km, double Rm) :PropensityFunction(REACTION_TYPE,1),s1(s1),s2(s2),s3(s3),Km(Km),Rm(Rm) {}
	uint s1;
	uint s2;
	uint s3;
	double Km;
	double Rm;

	void changeVolume(double volumeMultiplier) {}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		double u5 = double(speciesCounts[s1])/(Km + double(speciesCounts[s1]));
		return u5*Rm*(double(speciesCounts[s2]) - double(speciesCounts[s3]));
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderProductSubstrateDependent3Species needs one species dependencies, had", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "FirstOrderProductSubstrateDependent3Species needs one species dependencies, had", sd2.len);
		utuple sd3 = getSpecificDependencies(reactionIndex, D, 3);
		if (sd3.len != 1) throw InvalidArgException("D", "FirstOrderProductSubstrateDependent3Species needs one species dependencies, had", sd3.len);

		// Find the rate constant.
		if (k.len < 2) throw InvalidArgException("k", "FirstOrderProductSubstrateDependent3Species needs two parameters, had", k.len);

		return new FirstOrderProductSubstrateDependent3Species(sd1[0], sd2[0], sd3[0], k[0], k[1]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"item", "1/second", NULL};	
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderProductSubstrateDependent3Species", "(x1 / (k1 + x1)) * k2 * (x2 - x3)", unitsForConstants, &create);
	}
};

class FirstOrderU6 : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4007;

	FirstOrderU6(uint s1, double a1, double a2, double sc, double k1, double k2) :PropensityFunction(REACTION_TYPE,1),s1(s1),a1(a1),a2(a2),sc(sc),k1(k1),k2(k2) {}
	uint s1;
	double a1;
	double a2;
	double sc;
	double k1;
	double k2;

	void changeVolume(double volumeMultiplier) {}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		double u6 = 1/(1 + a1*exp(a2*(sc - double(speciesCounts[s1]))));
		return u6*k1*k2*double(speciesCounts[s1]);
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderU6 needs one species dependencies, had", sd1.len);
		
		// Find the rate constant.
		if (k.len < 5) throw InvalidArgException("k", "FirstOrderU6 needs five parameters, had", k.len);

		return new FirstOrderU6(sd1[0], k[0], k[1], k[2], k[3], k[4]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"", "1/item", "item", "", "1/second", NULL};
		const char* expressions[] = {"(1 / (1 + k1 * exp(k2 * (k3 - x1)))) * k4 * k5 * x1", "k4 * k5 * x1 * (1 / (1 + k1 * exp(k2 * (k3 - x1))))", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderU6", expressions, unitsForConstants, &create);
	}
};

class FirstOrderU7 : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4008;

	FirstOrderU7(uint s1, double a1, double a2, double sc, double k1):PropensityFunction(REACTION_TYPE,1),s1(s1),a1(a1),a2(a2),sc(sc),k1(k1) {}
	uint s1;
	double a1;
	double a2;
	double sc;
	double k1;

	void changeVolume(double volumeMultiplier) {}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		double u7 = 1/(1 + a1*exp(a2*(sc - double(speciesCounts[s1]))));
		return u7*k1*double(speciesCounts[s1]);
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderU7 needs one species dependencies, had", sd1.len);
		
		// Find the rate constant.
		if (k.len < 4) throw InvalidArgException("k", "FirstOrderU7 needs four parameters, had", k.len);

		return new FirstOrderU7(sd1[0], k[0], k[1], k[2], k[3]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"", "1/item", "item", "1/second", NULL};
		const char* expressions[] = {"(1 / (1 + k1 * exp(k2 * (k3 - x1)))) * k4 * x1", "k4 * x1 * (1 / (1 + k1 * exp(k2 * (k3 - x1))))", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderU7", expressions, unitsForConstants, &create);
	}
};

class SecondOrdercalc : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4009;

	SecondOrdercalc(uint s1, uint s2, double k1, double k2):PropensityFunction(REACTION_TYPE,2),s1(s1),s2(s2),k1(k1),k2(k2) {}
	uint s1;
	uint s2;
	double k1;
	double k2;

	void changeVolume(double volumeMultiplier) {k1/=volumeMultiplier;}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		return k1*double(speciesCounts[s1])*(k2 - double(speciesCounts[s2]));
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
                utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
                if (sd1.len != 1) throw InvalidArgException("D", "SecondOrdercalc needs one species dependencies, had", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "SecondOrdercalc needs one species dependencies, had", sd2.len);

		// Find the rate constant.
		if (k.len < 2) throw InvalidArgException("k", "SecondOrdercalc needs two parameters, had", k.len);

		return new SecondOrdercalc(sd1[0], sd2[0], k[0], k[1]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"1/(item*second)", "item", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "SecondOrdercalc", "k1 * x1 * (k2 - x2)", unitsForConstants, &create);
	}
};

class FourthOrdercalm : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4010;

	FourthOrdercalm(uint s1, uint s2, uint s3, double k1, double k2) :PropensityFunction(REACTION_TYPE,4),s1(s1),s2(s2),s3(s3),k1(k1),k2(k2) {}
	uint s1;
	uint s2;
	uint s3;
	double k1;
	double k2;

	void changeVolume(double volumeMultiplier) {k1/=(volumeMultiplier*volumeMultiplier*volumeMultiplier);}

	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		int s11 = speciesCounts[s1];
		int s22 = speciesCounts[s2];
		int s33 = speciesCounts[s3];		
		return k1 * double(s11*(s11-1)*(s11-2)) * (k2 - s22 - s33);
	}

#ifdef OPT_AVX
        avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
        {
                return naiveCalculateAvx(this, time, speciesCounts, numberSpecies);
        }
#endif

	static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
	{
		// Find the species dependency.
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FourthOrdercalm had invalid number of s1 dependencies", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "FourthOrdercalm had invalid number of s2 dependencies", sd2.len);
		utuple sd3 = getSpecificDependencies(reactionIndex, D, 3);
		if (sd3.len != 1) throw InvalidArgException("D", "FourthOrdercalm had invalid number of s3 dependencies", sd3.len);

		// Find the rate constants.
		if (k.len < 2) throw InvalidArgException("k", "FourthOrdercalm needs two parameters, had", k.len);

		return new FourthOrdercalm(sd1[0], sd2[0], sd3[0], k[0], k[1]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"1/(second*item*item*item)", "item", NULL};
		const char* expressions[] = {"k1 * x1^3 * (k2 - x2 - x3)", "k1 * x1 * x1 * x1 * (k2 - x2 - x3)", "k1 * x1 * (x1 - 1) * (x1 - 2) * (k2 - x2 - x3)", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FourthOrdercalm", expressions, unitsForConstants, &create);
	}
};

class FirstOrderMichaelisMentenRho : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4011;

	FirstOrderMichaelisMentenRho(uint s1, uint s2, uint s3, uint s4, double k1, double k2, double k3, double N, double lda) :PropensityFunction(REACTION_TYPE,1),s1(s1),s2(s2),s3(s3),s4(s4),k1(k1),k2(k2),k3(k3),N(N),lda(lda) {}
	uint s1;
	uint s2;
	uint s3;
	uint s4;
	double k1;
	double k2;
	double k3;
	double N;
	double lda;

	void changeVolume(double volumeMultiplier) {}
	double calculate (const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		double x = (1/double(speciesCounts[s1]))*((3*pow(double(speciesCounts[s2]),2))/(pow(k1,2) + pow(double(speciesCounts[s2]),2)))*(1 - ((0.8*pow(double(speciesCounts[s3]),2))/(pow(k2,2) + pow(double(speciesCounts[s3]),2))));
		double L0 = pow(10, -N/2);
		double xN1 = pow(x,N+1) - 1;
		double ldaxN1 = pow(lda*x,N+1) - 1;
		double multterm = (ldaxN1*(x-1))/(((lda*x)-1)*xN1);
		double rho = 1/(1+(L0*multterm));
		return k3*rho*double(speciesCounts[s4]);
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRho needs one species dependencies, had", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRho needs one species dependencies, had", sd2.len);
		utuple sd3 = getSpecificDependencies(reactionIndex, D, 3);
		if (sd3.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRho needs one species dependencies, had", sd3.len);
		utuple sd4 = getSpecificDependencies(reactionIndex, D, 4);
		if (sd4.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRho needs one species dependencies, had", sd4.len);

		// Find the rate constant.
		if (k.len < 5) throw InvalidArgException("k", "FirstOrderMichaelisMentenRho needs five parameters, had", k.len);

		return new FirstOrderMichaelisMentenRho(sd1[0], sd2[0], sd3[0], sd4[0], k[0],k[1], k[2], k[3], k[4]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"item", "item", "1/second", "", "", NULL};		
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderMichaelisMentenRho", "x4 * k3 * (1 / (1 + 10^(-k5 / 2) * ((k4 * 3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2))^(k5 + 1) - 1) * (3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2) - 1) / ((k4 * 3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2) - 1) * ((3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2))^(k5 + 1) - 1))))", unitsForConstants, &create);
	}
};

class FirstOrderMichaelisMentenRhoMinus : public lm::me::PropensityFunction
{
public:
        static const uint REACTION_TYPE = 4012;

        FirstOrderMichaelisMentenRhoMinus(uint s1, uint s2, uint s3, uint s4, double k1, double k2, double k3, double N, double lda) :PropensityFunction(REACTION_TYPE,1),s1(s1),s2(s2),s3(s3),s4(s4),k1(k1),k2(k2),k3(k3),N(N),lda(lda) {}
        uint s1;
        uint s2;
        uint s3;
        uint s4;
        double k1;
        double k2;
        double k3;
        double N;
        double lda;

        void changeVolume(double volumeMultiplier) {}
        double calculate (const double time, const int* speciesCounts, const uint numberSpecies) const
        {
                double x = (1/double(speciesCounts[s1]))*((3*pow(double(speciesCounts[s2]),2))/(pow(k1,2) + pow(double(speciesCounts[s2]),2)))*(1 - ((0.8*pow(double(speciesCounts[s3]),2))/(pow(k2,2) + pow(double(speciesCounts[s3]),2))));
                double L0 = pow(10, -N/2);
                double xN1 = pow(x,N+1) - 1;
                double ldaxN1 = pow(lda*x,N+1) - 1;
                double multterm = (ldaxN1*(x-1))/(((lda*x)-1)*xN1);
                double rho = 1/(1+(L0*multterm));
                return k3*(1-rho)*double(speciesCounts[s4]);
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
                utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
                if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRhoMinus needs one species dependencies, had", sd1.len);
                utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
                if (sd2.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRhoMinus needs one species dependencies, had", sd2.len);
                utuple sd3 = getSpecificDependencies(reactionIndex, D, 3);
                if (sd3.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRhoMinus needs one species dependencies, had", sd3.len);
                utuple sd4 = getSpecificDependencies(reactionIndex, D, 4);
                if (sd4.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenRhoMinus needs one species dependencies, had", sd4.len);

                // Find the rate constant.
                if (k.len < 5) throw InvalidArgException("k", "FirstOrderMichaelisMentenRhoMinus needs five parameters, had", k.len);

                return new FirstOrderMichaelisMentenRhoMinus(sd1[0], sd2[0], sd3[0], sd4[0], k[0],k[1], k[2], k[3], k[4]);
        }

        static lm::me::PropensityFunctionDefinition registerFunction()
        {
		const char* unitsForConstants[] = {"item", "item", "1/second", "", "", NULL};
                return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderMichaelisMentenRhoMinus", "x4 * k4 * (1 - 1 / (1 + 10^(-k3 / 2) * ((k5 * 3 * (1 - 0.8 * x1^2 / (k1^2 + x1^2)) * (1 / x3) * x2^2 / (k2^2 + x2^2))^(k3 + 1) - 1) * (3 * (1 - 0.8 * x1^2 / (k1^2 + x1^2)) * (1 / x3) * x2^2 / (k2^2 + x2^2) - 1) / ((k5 * 3 * (1 - 0.8 * x1^2 / (k1^2 + x1^2)) * (1 / x3) * x2^2 / (k2^2 + x2^2) - 1) * ((3 * (1 - 0.8 * x1^2 / (k1^2 + x1^2)) * (1 / x3) * x2^2 / (k2^2 + x2^2))^(k3 + 1) - 1))))", unitsForConstants, &create);
        }
};

class FirstOrderMichaelisMentenU2 : public lm::me::PropensityFunction
{
public:
	static const uint REACTION_TYPE = 4013;

	FirstOrderMichaelisMentenU2(uint s1, uint s2, uint s3, uint s4, uint s5, uint s6, double k1, double k2, double N, double lda, double Rm, double Km, double S) :PropensityFunction(REACTION_TYPE,1),s1(s1),s2(s2),s3(s3),s4(s4),s5(s5),s6(s6),k1(k1),k2(k2),N(N),lda(lda),Rm(Rm),Km(Km),S(S) {}
	uint s1;
	uint s2;
	uint s3;
	uint s4;
	uint s5;
	uint s6;
	double k1;
	double k2;
	double N;
	double lda;
	double Rm;
	double Km;
	double S;

	void changeVolume(double volumeMultiplier) {}
	double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
	{
		double x = (1/double(speciesCounts[s1]))*((3*pow(double(speciesCounts[s2]),2))/(pow(k1,2) + pow(double(speciesCounts[s2]),2)))*(1 - ((0.8*pow(double(speciesCounts[s3]),2))/(pow(k2,2) + pow(double(speciesCounts[s3]),2))));
		double L0 = pow(10, -N/2);
		double xN1 = pow(x,N+1) - 1;
		double ldaxN1 = pow(lda*x,N+1) - 1;
		double theta = 1/((xN1/(x-1)) + (L0*(ldaxN1/((lda*x)-1))));
		double h = double(speciesCounts[s4])/(double(speciesCounts[s4]) + double(speciesCounts[s6]));
		double u2 = h*theta;
		return u2*((S*Rm*double(speciesCounts[s5]))/(Km + double(speciesCounts[s5])));
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
		utuple sd1 = getSpecificDependencies(reactionIndex, D, 1);
		if (sd1.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU2 needs one species dependencies, had", sd1.len);
		utuple sd2 = getSpecificDependencies(reactionIndex, D, 2);
		if (sd2.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU2 needs one species dependencies, had", sd2.len);
		utuple sd3 = getSpecificDependencies(reactionIndex, D, 3);
		if (sd3.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU2 needs one species dependencies, had", sd3.len);
		utuple sd4 = getSpecificDependencies(reactionIndex, D, 4);
		if (sd4.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU2 needs one species dependencies, had", sd4.len);
		utuple sd5 = getSpecificDependencies(reactionIndex, D, 5);
		if (sd5.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU2 needs one species dependencies, had", sd5.len);
		utuple sd6 = getSpecificDependencies(reactionIndex, D, 6);
		if (sd6.len != 1) throw InvalidArgException("D", "FirstOrderMichaelisMentenU2 needs one species dependencies, had", sd6.len);

		// Find the rate constant.
		if (k.len < 7) throw InvalidArgException("k", "FirstOrderMichaelisMentenU2 needs seven parameters, had", k.len);

		return new FirstOrderMichaelisMentenU2(sd1[0], sd2[0], sd3[0], sd4[0], sd5[0], sd6[0], k[0], k[1], k[2], k[3], k[4], k[5], k[6]);
	}

	static lm::me::PropensityFunctionDefinition registerFunction()
	{
		const char* unitsForConstants[] = {"item", "item", "", "", "item/second", "item", "", NULL};
		return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderMichaelisMentenU2", "k7 * k5 * x5 * (x4 / (x4 + x6)) * (1 / (((3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2))^(k3 + 1) - 1) / (3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2) - 1) + (10^(-k3/2)) * (((k4 * 3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2))^(k3 + 1) - 1) / (k4 * 3 * (1 - 0.8 * x3^2 / (k2^2 + x3^2)) * (1 / x1) * x2^2 / (k1^2 + x2^2) - 1)))) / (k6 + x5)", unitsForConstants, &create);
	}
};

list<lm::me::PropensityFunctionDefinition> EnzymePropensityFunctions::getPropensityFunctionDefinitions()
{
    list<lm::me::PropensityFunctionDefinition> defs;
    defs.push_back(FirstOrderMichaelisMenten::registerFunction());
    defs.push_back(FirstOrderMichaelisMentenU1::registerFunction());
    defs.push_back(FirstOrderDoubleMichaelisMenten::registerFunction());
    defs.push_back(FirstOrderProductSubstrateDependent2Species::registerFunction());
    defs.push_back(FirstOrderProductSubstrateDependent3Species::registerFunction());
    defs.push_back(FirstOrderU6::registerFunction());
    defs.push_back(FirstOrderU7::registerFunction());
    defs.push_back(SecondOrdercalc::registerFunction());
    defs.push_back(FourthOrdercalm::registerFunction());
    defs.push_back(FirstOrderMichaelisMentenRho::registerFunction());
    defs.push_back(FirstOrderMichaelisMentenRhoMinus::registerFunction());
    defs.push_back(FirstOrderMichaelisMentenU2::registerFunction());
    return defs;
}

}
}
