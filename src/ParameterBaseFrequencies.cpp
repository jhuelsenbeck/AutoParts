#include <cmath>
#include <iostream>
#include "ParameterBaseFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"

#define MIN_FREQ 1e-6



ParameterBaseFrequencies::ParameterBaseFrequencies(const ParameterBaseFrequencies& parm) : Parameter(parm) {

    std::cout << "ParameterBaseFrequencies copy constructor" << std::endl;

    freqs[0] = parm.freqs[0];
    freqs[1] = parm.freqs[1];
}

ParameterBaseFrequencies::ParameterBaseFrequencies(Model* m, UserSettings* s) : Parameter(m, s, "Base Frequencies") {

    updateModifiesEigens = true;
    freqs[0].resize(4);
    freqs[1].resize(4);
    std::vector<double> alpha(4, 1.0);
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    bool err = false;
    do
        {
        err = Probability::Dirichlet::rv(&rng, alpha, freqs[0]);
        } while (anyLessThanMin(freqs[0], MIN_FREQ) == true || err == true);
    freqs[1] = freqs[0];
}

void ParameterBaseFrequencies::accept(void) {

    freqs[1] = freqs[0];
}

std::string ParameterBaseFrequencies::getHeader(void) {

    return "";
}

std::vector<double> ParameterBaseFrequencies::getValue(void) {

    std::vector<double> vals;
    for (int i=0; i<freqs[0].size(); i++)
        vals.push_back(freqs[0][i]);
    return vals;
}

std::string ParameterBaseFrequencies::getValuesAsString(int precision) {

    std::string str = "(";
    for (int i=0; i<freqs[0].size(); i++)
        {
        str += std::to_string(freqs[0][i]);
        if (i+1 != freqs[0].size())
            str += ", ";
        }
    str += ")";
    return str;
}

double ParameterBaseFrequencies::lnProbability(void) {

    return log(6.0);
}

Parameter* ParameterBaseFrequencies::newRandomInstance(void) {

    return new ParameterBaseFrequencies(modelPtr, userSettings);
}

void ParameterBaseFrequencies::print(void) {

    std::cout << getValuesAsString(3) << std::endl;
}

void ParameterBaseFrequencies::reject(void) {

    freqs[0] = freqs[1];
}

std::string ParameterBaseFrequencies::type(void) {

    return "base frequencies";
}

double ParameterBaseFrequencies::update(void) {


    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    if(rng.uniformRv() < 0.5)
        {
        double alpha0 = UpdateInfo::updateInfo().attempt(UpdateType::BASE_FREQUENCY_DIRICHLET);
        std::vector<double> oldFreqs(4);
        std::vector<double> alphaForward(4);
        for (int i=0; i<4; i++)
            {
            oldFreqs[i] = freqs[0][i];
            alphaForward[i] = freqs[0][i] * alpha0;
            }
            
        std::vector<double> newFreqs(4);
        bool err = Probability::Dirichlet::rv(&rng, alphaForward, newFreqs);
        if (err == true || anyLessThanMin(newFreqs, MIN_FREQ) == true)
            return rejectProposal();     // outside the support (see Parameter.hpp)
        
        std::vector<double> alphaReverse(4);
        for (int i=0; i<4; i++)
            {
            alphaReverse[i] = newFreqs[i] * alpha0;
            freqs[0][i] = newFreqs[i];
            }
            
        double lnHastingsRatio = Probability::Dirichlet::lnPdf(alphaReverse, oldFreqs) - Probability::Dirichlet::lnPdf(alphaForward, newFreqs);
        return lnHastingsRatio;
    }
    else 
        {
        double lnHastingsRatio = 0.0;

        int randomIndex = (int)(rng.uniformRv() * 4);
        double alpha0 = UpdateInfo::updateInfo().attempt(UpdateType::BASE_FREQUENCY_BETA);

        double currentValue = freqs[0][randomIndex];
        std::vector<double> newFreqs(4);

        double a = alpha0 + 1.0;
        double b = alpha0 / currentValue - a + 2.0;
        {
            newFreqs = freqs[0];
            double newValue = Probability::Beta::rv(&rng, a, b);
            newFreqs[randomIndex] = newValue;

            double scaling = (1.0 - newValue) / (1.0 - currentValue);

            double sum = 0.0;
            for(int i = 0; i < 4; i++)
                {
                if(i != randomIndex)
                    {
                    newFreqs[i] *= scaling;
                    }

                sum += newFreqs[i];
                }

            for(int i = 0; i < 4; i++)
                {
                newFreqs[i] /= sum;
                }
            
            double newB = alpha0 / newFreqs[randomIndex] - a + 2.0;
            lnHastingsRatio = Probability::Beta::lnPdf(a, newB, currentValue) - Probability::Beta::lnPdf(a, b, newValue);
            lnHastingsRatio += 2 * std::log(scaling) - 3 * std::log(sum);
            }
        if (anyLessThanMin(newFreqs, MIN_FREQ) == true)
            return rejectProposal();     // outside the support (see Parameter.hpp)

        freqs[0] = newFreqs;

        return lnHastingsRatio;
        }
}
