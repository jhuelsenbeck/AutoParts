#include <cmath>
#include <iostream>
#include "ParameterTreeLength.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"

#define MIN_LENGTH 1e-6
#define MAX_LENGTH 50.0


ParameterTreeLength::ParameterTreeLength(const ParameterTreeLength& parm) : Parameter(parm) {

    std::cout << "ParameterTreeLength copy constructor" << std::endl;

    alphaT = parm.alphaT;
    betaT = parm.betaT;
    length[0] = parm.length[0];
    length[1] = parm.length[1];
}

ParameterTreeLength::ParameterTreeLength(Model* m, UserSettings* s, double a, double b) : Parameter(m, s, "Tree Length") {

    updateModifiesEigens = false;
    alphaT = a;
    betaT = b;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    do
        {
        length[0] = Probability::Gamma::rv(&rng, alphaT, betaT);
        } while (length[0] < MIN_LENGTH || length[0] > MAX_LENGTH);
    length[1] = length[0];
}

void ParameterTreeLength::accept(void) {

    length[1] = length[0];
}

std::string ParameterTreeLength::getHeader(void) {

    return "";
}

std::vector<double> ParameterTreeLength::getValue(void) {

    std::vector<double> vals;
    vals.push_back(length[0]);
    return vals;
}

std::string ParameterTreeLength::getValuesAsString(int precision) {

    std::string str = "";
    str += std::to_string(length[0]);
    return str;
}

double ParameterTreeLength::lnProbability(void) {
    double c = Probability::Gamma::cdf(alphaT, betaT, MAX_LENGTH) - Probability::Gamma::cdf(alphaT, betaT, MIN_LENGTH);
    return (Probability::Gamma::lnPdf(alphaT, betaT, length[0]) - log(c));
}

Parameter* ParameterTreeLength::newRandomInstance(void) {

    return new ParameterTreeLength(modelPtr, userSettings, alphaT, betaT);
}

void ParameterTreeLength::print(void) {

    std::cout << getValuesAsString(3) << std::endl;
}

void ParameterTreeLength::reject(void) {

    length[0] = length[1];
}

std::string ParameterTreeLength::type(void) {

    return "tree length";
}

double ParameterTreeLength::update(void) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
        double tuning = UpdateInfo::updateInfo().attempt(UpdateType::TREE_LENGTH);
    double randomFactor = exp( (rng.uniformRv()-0.5)*tuning );
    double newVal = length[0] * randomFactor;
    if (newVal < MIN_LENGTH || newVal > MAX_LENGTH)
        return rejectProposal();     // outside the support (see Parameter.hpp)
    length[0] = newVal;
    return log(randomFactor);
}
