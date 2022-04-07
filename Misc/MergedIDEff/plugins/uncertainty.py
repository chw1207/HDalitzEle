import numpy as np

# This is the script to estimate the uncertainty of efficiency.
# Several methods are listed below. (materials i can find for now)
# No idea which method should be used for now........ 
# Fortunately, they do not give much difference when statistic is enough.

# Some references:
#*  1. Treatment of Errors in Efficiency Calculations (Brookhaven's paper)
#*  2. Efficiency measurement: a Bayesian approach (NYU's paper)
#*  3. TEfficiency document: https://root.cern.ch/doc/master/classTEfficiency.html
#*  4. https://indico.cern.ch/event/66256/contributions/2071577/attachments/1017176/1447814/EfficiencyErrors.pdf 

# The most safe way? -> use ROOT package
#*  1. TEfficiency and TGraphAsymmErrors.BayesDivide()
#*  2. TEfficiency.SetStatisticOption(ROOT.TEfficiency.kBUniform)
#*  3. TEfficiency.SetConfidenceLevel(0.683)
#*  4. Two methods give the same value of efficiency and uncertainty
#*  5. HOWEVER, both of them only work for 1D efficiency -> how about 2D efficiency? 

# Uncertainty of the Efficiency
#* Some notes:
#*  1. The process is assumed to be the Bernoulli experiment.
#*  2. The efficiency is then the probability for an experiment to give a 
#*     positive outcome.(eff. = # of successes / # of events obtained by 
#*     an offline reconstruction)
#*  3. Then the Binomial error(eff) = sqrt(eff * (1.-eff)/ # of events 
#*     obtained by an offline reconstruction)

# Most widely used way
# Leads to unphysical reselts when numerator = 0 and numerator = denominator (eff = numerator / denominator)
# CMS-muon-pog: https://gitlab.cern.ch/cms-muonPOG/spark_tnp/-/blob/master/prepare.py#L25-29
def Binomial_Error(eff, N):
    error = np.sqrt(eff*(1.-eff)/N)
    return error


# Another approch found at refernce[1, 2]
# Solves the problem of binomial error leading to unphysical results
def Bayesian_Error(num, den):
    k, n = num, den
    error = np.sqrt(((k+1)*(k+2))/((n+2)*(n+3)) - np.power(k+1, 2)/np.power(n+2, 2))
    return error


# Error Propagation for Scale Factor (SF = Eff_A / Eff_B)
#* Some notes:
#*  1. When the quantity which depends on other measurements, the error
#*     of this quantity comes from other measurements.
#*  2. In other words, the uncertainties of other measurements propagate to
#*     the uncertainty of this quantity.
#*  3. For example, the uncertainty of the scale factor comes from 
#*     the uncertainty of the effs. of data and MC.
def Error_Propagation(f, sigmaA, A, sigmaB, B):
    error = f * np.sqrt((sigmaA/A)*(sigmaA/A)+(sigmaB/B)*(sigmaB/B))
    return error
