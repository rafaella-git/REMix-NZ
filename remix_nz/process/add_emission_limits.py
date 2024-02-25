import pandas as pd

def add_emission_limit(m):
    # Add net zero restriction for 2050
    accounting_emissionLimit = pd.DataFrame(
        index=pd.MultiIndex.from_product([["global"], [2050], ["CO2Emission"]])
    )
    accounting_emissionLimit["useUpper"] = 1  # minimization of system costs
    accounting_emissionLimit["upperValue"] = 0  # minimization of system costs
    m["Base"].parameter.add(accounting_emissionLimit, "accounting_indicatorbounds")

def add_emission_budget(m):
    # Add cumulativ emission budget for all years
    accounting_emissionBudget = pd.DataFrame(
        index=pd.MultiIndex.from_product([["global"], ["horizon"], ["CO2Emission"]])
    )
    accounting_emissionBudget["integral"] = 1
    accounting_emissionBudget["endyear"] = 25 # length of the last year to run (2050)
    accounting_emissionBudget["useUpper"] = 1
    accounting_emissionBudget["upperValue"] = 45000  # 45 Gt CO2
    m["Base"].parameter.add(accounting_emissionBudget, "accounting_indicatorbounds")
