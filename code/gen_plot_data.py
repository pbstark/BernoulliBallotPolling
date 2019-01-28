import sys
import math
import scipy as sp
import scipy.stats
from scipy.stats import norm
from scipy.stats import poisson, binom, hypergeom
from permute.utils import binom_conf_interval, hypergeom_conf_interval
from numpy.random import choice
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sprt
from joblib import Parallel, delayed
import multiprocessing
import time

num_cores = multiprocessing.cpu_count()

def run_bpa_trial(prop_winner, Ntot, alpha, trial_num):
    if not trial_num%10000:
        print trial_num

    test = 1
    c_samp = 0
    sample_prob = (prop_winner/float(Ntot))
    votes = sp.stats.binom.rvs(1, sample_prob, size=Ntot)

    ctr = 0
    for vote in votes:
        if test >= 1/alpha:
            break
        if vote:
            test = test*sample_prob/.5
        else:
            test = test*(1 - sample_prob)/.5

        c_samp += 1

    return c_samp

def get_bpa_sample_size(prop_winner, Ntot, alpha):
    n_trials = 10**7
    trials = []

    trials = Parallel(n_jobs=num_cores)(delayed(run_bpa_trial)(prop_winner, Ntot, alpha, i) \
            for i in range(n_trials))

    sorted(trials)
    quants = {}
    for quant in [25, 50, 75, 90, 99]:
        quants[quant] = np.percentile(trials, quant)/float(Ntot)

    return quants 

def compute(margin, N_wl, pi, alpha, n): #start, end):

    ret_val = 0
    p1 = 0.5 + margin/2
    q_alpha = hypergeom.ppf(1-alpha, N_wl, N_wl/2, n) # upper alpha quantile of null distribution
    prob_select_N = binom.pmf(n, N_wl, pi) # probability of selecting n out of N_wl at sampling rate pi
    pvalue_nw = hypergeom.sf(q_alpha, N_wl, N_wl*p1, n) # probability of alternative distr falling above q_alpha
    return prob_select_N*pvalue_nw

def compute_unconditional_power(margin, N_wl, pi, alpha):
    '''
    Compute unconditional power of the test.

    margin = vote margin (votes for w / votes for w or l) in the population
    N_wl = the total number of ballots for either the winner or loser in the population,
    pop = total population size,
    pi = the sampling probability,
    alpha = the type I error rate
    '''
    unlikely_draw_lower = binom.ppf(0.005, N_wl, pi)
    unlikely_draw_upper = binom.ppf(0.995, N_wl, pi)
    power_sum = 0

    powers = Parallel(n_jobs=num_cores)(delayed(compute)(margin, N_wl, pi, alpha, n) \
            for n in range(int(unlikely_draw_lower), int(unlikely_draw_upper)))
    return sum(powers)

def get_sample_for_power(margin, Ntot, alpha, power, base_rate):
    pi = 0
    while True:
        pi += base_rate 
        if pi > 1:
            return 1
        x = compute_unconditional_power(margin, Ntot, pi, alpha)
        if x >= power:
            return pi

def get_bbp_sample_size(prop_winner, Ntot, alpha):
    quants = {}
    for quant in [25, 50, 75, 90, 99]:
        margin = (prop_winner - (Ntot - prop_winner))/float(Ntot)
        print quant
        quants[quant] = get_sample_for_power(margin, Ntot, alpha, quant/100.0, 1/float(Ntot))

    return quants#, ASN

def main():
    prop_ws = []

    for i in range(1, 50):
        prop_ws.append((100-i)/100.0)

    prop_ws.reverse()

    data_file =  open('data.csv', 'w')
    data_file.write('alpha,total,prop_w,bpa_25,bpa_50,bpa_75,bpa_90,bpa_99,bbp_25,bbp_50,bbp_75,bbp_90,bbp_99\n')
    idx = 0

    for alpha in [.1]:#, .05, .01]:#.05, .1, .01]:
        Ntot = 10000
        bbp = []
        bpa = []
        for prop_w in prop_ws:
            print prop_w
            start_time = time.time()
            prop_winner = Ntot*prop_w
            bpa_ss = get_bpa_sample_size(prop_winner, Ntot, alpha)
            bbp_ss = get_bbp_sample_size(prop_winner, Ntot, alpha)

            print('Time: {}s'.format(time.time() - start_time))

            #save our data for later
            out_line = '{},{},{},'.format(alpha, Ntot, prop_w)
            bpa_string = '{},{},{},{},{},'.format(bpa_ss[25],bpa_ss[50],bpa_ss[75],bpa_ss[90],bpa_ss[99])
            bbp_string = '{},{},{},{},{}'.format(bbp_ss[25],bbp_ss[50],bbp_ss[75],bbp_ss[90],bbp_ss[99])
            out_line += bpa_string + bbp_string + '\n'
            data_file.write(out_line)
            data_file.flush()

            bpa.append(bpa_ss)
            bbp.append(bbp_ss)

        print('Computed 10k')

	start_time = time.time()
        Ntot = 1000000
        bbp_1m = []
        bpa_1m = []
        for prop_w in prop_ws:
            print '1M prop_w:',  prop_w
            prop_winner = Ntot*prop_w
            bpa_ss = get_bpa_sample_size(prop_winner, Ntot, alpha)
            bbp_ss = get_bbp_sample_size(prop_winner, Ntot, alpha)

            #save our data for later
            out_line = '{},{},{},'.format(alpha, Ntot, prop_w)
            bpa_string = '{},{},{},{},{},'.format(bpa_ss[25],bpa_ss[50],bpa_ss[75],bpa_ss[90],bpa_ss[99])
            bbp_string = '{},{},{},{},{}'.format(bbp_ss[25],bbp_ss[50],bbp_ss[75],bbp_ss[90],bbp_ss[99])
            out_line += bpa_string + bbp_string + '\n'
            data_file.write(out_line)
            data_file.flush()
            bpa_1m.append(bpa_ss)
            bbp_1m.append(bbp_ss)

        print('Computed 1M in {} seconds'.format(time.time() - start_time))

            #idx += 1
    data_file.close()
if __name__== "__main__":
      main()
