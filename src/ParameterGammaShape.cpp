#include <cmath>
#include <iostream>
#include "ParameterGammaShape.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"

#define MIN_SHAPE 1e-6
#define MAX_SHAPE 10.0



ParameterGammaShape::ParameterGammaShape(const ParameterGammaShape& parm) : Parameter(parm) {

    std::cout << "ParameterGammaShape copy constructor" << std::endl;

    lambda = parm.lambda;
    shape[0] = parm.shape[0];
    shape[1] = parm.shape[1];
    numGammaCategories = parm.numGammaCategories;
    gammaRates[0] = parm.gammaRates[0];
    gammaRates[1] = parm.gammaRates[1];
}

ParameterGammaShape::ParameterGammaShape(Model* m, UserSettings* s, double lam, int nc) : Parameter(m, s, "Gamma Shape") {

    updateModifiesEigens = false;
    lambda = lam;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    do
        {
        shape[0] = Probability::Exponential::rv(&rng, lambda);
        } while (shape[0] < MIN_SHAPE || shape[0] > MAX_SHAPE);
    shape[1] = shape[0];

    numGammaCategories = nc;
    gammaRates[0].resize(numGammaCategories);
    Probability::Gamma::discretization(gammaRates[0], shape[0], shape[0], numGammaCategories, false);
    gammaRates[1] = gammaRates[0];
}

void ParameterGammaShape::accept(void) {

    shape[1] = shape[0];
    gammaRates[1] = gammaRates[0];
}

std::string ParameterGammaShape::getHeader(void) {

    return "";
}

std::vector<double> ParameterGammaShape::getValue(void) {

    std::vector<double> vals;
    vals.push_back(shape[0]);
    return vals;
}

std::string ParameterGammaShape::getValuesAsString(int precision) {

    std::string str = "";
    str += std::to_string(shape[0]);
    return str;
}

double ParameterGammaShape::lnProbability(void) {

    double c = Probability::Exponential::cdf(lambda, MAX_SHAPE) - Probability::Exponential::cdf(lambda, MIN_SHAPE);
    return Probability::Exponential::lnPdf(lambda, shape[0]) - log(c);
}

Parameter* ParameterGammaShape::newRandomInstance(void) {

    return new ParameterGammaShape(modelPtr, userSettings, lambda, numGammaCategories);
}

void ParameterGammaShape::print(void) {

    std::cout << getValuesAsString(3) << std::endl;
}

void ParameterGammaShape::reject(void) {

    shape[0] = shape[1];
    gammaRates[0] = gammaRates[1];
}

std::string ParameterGammaShape::type(void) {

    return "gamma shape";
}

double ParameterGammaShape::update(void) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double tuning = UpdateInfo::updateInfo().attempt(UpdateType::SHAPE_SCALE);
    
    double randomFactor = exp( (rng.uniformRv()-0.5)*tuning );
    double newVal = shape[0] * randomFactor;
    if (newVal < MIN_SHAPE || newVal > MAX_SHAPE)
        return rejectProposal();     // outside the support (see Parameter.hpp)
    
    shape[0] = newVal;
    Probability::Gamma::discretization(gammaRates[0], shape[0], shape[0], numGammaCategories, false);

    return log(randomFactor);
}
