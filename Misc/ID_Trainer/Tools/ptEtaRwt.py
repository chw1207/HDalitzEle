import pandas as pd


# deal with pT-eta reweighting
# pt eta will be reweighted to the class (cand)
def df_pteta_rwt(
    Mdf, label,
    ptw=[10, 30, 40, 50, 200, 10000], etaw=[-1.5, -1.0, 1.0, 1.5], eta="", pt="",
    SumWeightCol="instwei", NewWeightCol="NewWt", cand="", Classes=[""]
) -> pd.Series:

    Mdf["rwt"] = 1.
    Mdf[NewWeightCol] = 1.

    ptwt = [1.0] * len(ptw)
    etawt = [1.0] * len(etaw)

    for k in range(len(etaw)):
        if k == len(etaw) - 1:
            continue
        for i in range(len(ptw)):
            if i == len(ptw) - 1:
                continue
            for target in Classes:
                if target != cand:
                    targetSum = Mdf.loc[(Mdf[pt] <= ptw[i+1]) & (Mdf[pt] > ptw[i]) & (Mdf[eta] <= etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == target), SumWeightCol].sum()
                    candSum = Mdf.loc[(Mdf[pt] <= ptw[i+1]) & (Mdf[pt] > ptw[i]) & (Mdf[eta] <= etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == cand), SumWeightCol].sum()
                    if candSum > 0 and targetSum > 0:
                        ptwt[i] = candSum / (targetSum)
                    else:
                        ptwt[i] = 0

                    Mdf.loc[(Mdf[pt] <= ptw[i+1]) & (Mdf[pt] > ptw[i]) & (Mdf[eta] <= etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == cand), "rwt"] = 1.0
                    Mdf.loc[(Mdf[pt] <= ptw[i+1]) & (Mdf[pt] > ptw[i]) & (Mdf[eta] <= etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == target), "rwt"] = ptwt[i]

    Mdf.loc[:, NewWeightCol] = Mdf.loc[:, "rwt"] * Mdf.loc[:, SumWeightCol]

    for justclass in Classes:
        Sum = Mdf.loc[Mdf[label] == justclass, NewWeightCol].sum()
        print("Number of events in {} after pt-eta reweighting = {:.2f}".format(justclass, Sum), flush=True)

    return Mdf[NewWeightCol]


# balanced reweighting
# unbalanced multi-classification: https://datascience.stackexchange.com/a/49067
def df_balance_rwt(Mdf, SumWeightCol="instwei", NewWeightCol="NewWt", Classes=[""]) -> pd.Series:
    Mdf[NewWeightCol] = 1
    sum_w, wei = [1.0] * len(Classes), [1.0] * len(Classes)
    for i, k in enumerate(Classes):
        sum_w[i] = Mdf[SumWeightCol][Mdf.Class == k].sum()
    minimum = min(sum_w)

    for i, k in enumerate(Classes):
        wei[i] = minimum / sum_w[i]
        Mdf.loc[Mdf.Class == k, "rwt"] = wei[i]
        Mdf.loc[:, NewWeightCol] = Mdf.loc[:, "rwt"] * Mdf.loc[:, SumWeightCol]
        print("Class = %s, n = %.2f, balanced weight = %.2e" %(k, sum_w[i], wei[i]), flush=True)

    return Mdf[NewWeightCol]