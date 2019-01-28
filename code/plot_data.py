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
import csv

num_cores = multiprocessing.cpu_count()

if len(sys.argv) != 2:
    print('usage: python plot_data.py [data_file].csv')
    sys.exit()

def main():

    data = {}

    for line in csv.DictReader(open(sys.argv[1])):
        alpha = float(line['alpha'])
        prop_w = float(line['prop_w'])
        total = int(line['total'])
        to_set = { 
            'bpa_75': float(line['bpa_75']),
            'bpa_90': float(line['bpa_90']),
            'bpa_99': float(line['bpa_99']),
            'bbp_75': float(line['bbp_75']),
            'bbp_90': float(line['bbp_90']),
            'bbp_99': float(line['bbp_99']),
            'bbp_seq_75': float(line['bbp_seq_75']),
            'bbp_seq_90': float(line['bbp_seq_90']),
            'bbp_seq_99': float(line['bbp_seq_99']),
            }


        if alpha in data:
            if prop_w in data[alpha]:
                data[alpha][prop_w][total] = to_set
            else:
                data[alpha][prop_w] = {total: to_set}
        else:
            data[alpha] = {prop_w:  {total: to_set}}

        if data[alpha][prop_w][total]['bbp_90'] == None:
            data[alpha][prop_w][total]['bbp_90'] = .01
    prop_ws = []

    for i in range(1, 50):
        prop_ws.append((100-i)/100.0)

    prop_ws.reverse() 

    cols = ['Alpha {}'.format(col) for col in sorted(data.keys(), reverse=True)]
    rows = ['Quant {}'.format(row) for row in ['75', '90', '99']]
    fig, axes = plt.subplots(nrows=3, ncols=3)
       
    i = 1
    tots = [10000, 1000000]
    for quant in ['75', '90', '99']:
        #ax1 = plt.subplot(5, 3, 0)
        for alpha in data:
            bpa = []
            bbp = []
            bbp_seq = []
            bpa_1m = []
            bbp_1m = []
            bbp_seq_1m = []
            for prop_w in prop_ws: 
                bpa.append(data[alpha][prop_w][10000]['bpa_' + quant])
                bbp.append(data[alpha][prop_w][10000]['bbp_' + quant])
                bbp_seq.append(data[alpha][prop_w][10000]['bbp_seq_' + quant])
                bpa_1m.append(data[alpha][prop_w][1000000]['bpa_' + quant])
                bbp_1m.append(data[alpha][prop_w][1000000]['bbp_' + quant])
                bbp_seq_1m.append(data[alpha][prop_w][1000000]['bbp_seq_' + quant])
            
            ax1 = plt.subplot(3, 3, i)
            ax1.invert_xaxis()

            #plt.title('Quant: ' + str(quant) + ' Alpha: ' + str(alpha), fontsize=12)
            l1 = ax1.semilogy(prop_ws, bpa, label='BPA10k')
            l2 = ax1.semilogy(prop_ws, bbp, label='BBP10k')
            l3 = ax1.semilogy(prop_ws, bbp_seq, label='BBPSEQ10k')

            l3 = ax1.semilogy(prop_ws, bpa_1m, label='BPA1M')
            l4 = ax1.semilogy(prop_ws, bbp_1m, label='BBP1M')
            l4 = ax1.semilogy(prop_ws, bbp_seq_1m, label='BBPSEQ1M')

            if i%3 != 0:
                ax1.set_yticklabels([])
            else:
                ax1.yaxis.tick_right()

            if i%3 == 1:
                 ax1.set_ylabel(rows[i/3], rotation=90)
            if i < 7:
                ax1.set_xticklabels([])
            if i < 4:
                ax1.set_title(cols[i-1])
            if i == 8:
                ax1.set_xlabel('Margin', size='large')

            i += 1
            
            

    
    fig.text(0.04, 0.5, '% ballots', va='center', rotation='vertical', size='large')
    handles, labels = ax1.get_legend_handles_labels()
    art = []
    lgd = ax1.legend(handles, labels, loc=9, ncol=6, fontsize=10, bbox_to_anchor=(-.75, -.35))
    art.append(lgd)
    plt.savefig('quant_plot.png', additional_artists=art, bbox_inches="tight")



if __name__== "__main__":
      main()
