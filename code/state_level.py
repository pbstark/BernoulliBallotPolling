from __future__ import (absolute_import, division,
                                print_function, unicode_literals)
import math
import numpy as np
import scipy as sp
from scipy.stats import poisson, binom
from sprt import ballot_polling_sprt
import sys
import csv

def compute_power(Ntot, Vw, Vl, pi, alpha, reps=10**3):
    Vu = Ntot - Vw - Vl
    pop = np.array([1]*Vw + [0]*Vl + [np.nan]*Vu)
    
    power_sum = 0
    for i in range(reps):
        np.random.shuffle(pop)
        n = binom.rvs(Ntot, pi)
#        print(i, n)
        sam = pop[:n]
        res = ballot_polling_sprt(sample=sam, popsize=Ntot, 
                                  alpha=alpha, Vw=Vw, Vl=Vl)
        if res['pvalue']<=alpha:
            power_sum += 1
    return power_sum/reps

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 state_level.py [state_level data.csv] [sample probability] [risk limit]") 
        sys.exit()

    print("|{:>5} | {:>7} | {:5} | {:>5}| {:>7}| {:>7}|".format("State", "Total", "Margin",  "Power", "Workload", "BPA Workload"))
    print("|---:|---:|---:|---:|---:|---:|---:|")
    sample_prob = float(sys.argv[2]) 
    risk_limit = float(sys.argv[3])

    for item in csv.DictReader(open(sys.argv[1])):
        Ntot = int(item["total"].replace(",", ""))
        margin = (float(item["winner"]) - float(item["runner-up"]))/Ntot
        Vw = int(item["winner"])
        Vl = int(item["runner-up"])


        res = compute_power(Ntot, Vw, Vl, sample_prob, risk_limit, reps=10**3)

        p_w = float(item["winner"])/Ntot
        p_l = float(item["runner-up"])/Ntot
        s_w = p_w/(p_w + p_l)
        z_w = math.log(2*s_w)  
        z_l = math.log(2 - 2*s_w)
        
        bpa_asn =  (math.log(1.0/risk_limit) + z_w*0.5)/(p_w*z_w + p_l*z_l)

        print("|{:>5} | {:>7} | {:>2.4f} | {:>.2f} | {:>8} | {:>8} |"\
                .format(item["label"], Ntot, margin, res, int(Ntot*sample_prob), int(bpa_asn)))

main()
