""" import libraries """

import csv
import os
import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt
import seaborn as sns
import random

""" functions """

def bootrstrap(data):                   # generates a randomized dataset
    simul = []
    
    for n in range(N):
        temp_data = []
        for i in range(S):
            if perm_stages:
                temp_data.append(data[random.randint(0, N - 1)][random.randint(0, S - 1)])
            else:
                temp_data.append(data[random.randint(0, N - 1)][i])
        
        simul.append(temp_data)
    
    return np.array(simul)


def normalize(data):                    # returns fractional data
    norm_array = []
    dataset = data.transpose()
    
    for trace in dataset:
        norm = []
        total = np.sum(trace)
        
        for it in trace:
            norm.append(it / total)
        
        norm_array.append(norm)
    
    out = np.array(norm_array)
    
    return out.transpose()


def cvdiff(data):                       # generates the difference between absolute and fractional CoV
    absolute = [] 
    for i in range(S):
       absolute.append(sts.tstd(data[i]) / sts.tmean(data[i]))
      
    norm_data = normalize(data)
    relative = []
    for i in range(S):
       relative.append(sts.tstd(norm_data[i]) / sts.tmean(norm_data[i]))
    
    delta = []
    for i in range(S):
        delta.append(absolute[i] - relative[i])
    
    return delta


def significance(simuls):               # returns the minimal significance level based of sig
    delt = simuls.transpose()
    delt = np.sort(delt)
    sig_i1 = int(bootstraps * sig1)
    if more_sig:
        sig_i2 = int(bootstraps * sig2)
        sig_i3 = int(bootstraps * sig3)
    
    sigs_1 = []
    if more_sig:
        sigs_2 = []
        sigs_3 = []
    for s in range(S):
        sigs_1.append([delt[s][len(delt[s]) - sig_i1 - 1]])
        if more_sig:
            sigs_2.append([delt[s][len(delt[s]) - sig_i2 - 1]])
            sigs_3.append([delt[s][len(delt[s]) - sig_i3 - 1]])
    
    if not more_sig:
        return sigs_1
    else:
        return sigs_1, sigs_2, sigs_3


def pval(emp, simuls):                  # returns empirical p-values for CVabs-CVfrac differences
    pvals = []
    
    for i in range(S):
        cnt = 0
        
        for j in range(bootstraps):
            if simuls[j][i] >= emp[i]:
                cnt += 1
        
        pvals.append(cnt / bootstraps)
    
    return pvals


""" configs """

source_file = 'thaliana_growth'         # file name used to take the raw data from (up to '.csv')
bootstraps = 10_000                     # number of bootstraps to run
more_sig = False                        # plot further significance points
sig1 = 0.05                             # significance threshold 1
sig2 = 0.01                             # significance threshold 2 
sig3 = 0.001                            # significance threshold 3 
save_boots = False                      # save all bootsraps in Excel
save_sigs = False                       # save empirical significances in Excel
perm_stages = False                     # bootsrtap individuals AND stages

""" start """

directory = os.getcwd()

read = []
with open(os.path.join(directory, source_file + '.csv'), newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        tmp = []
        for item in row:
            tmp.append(float(item))
        read.append(tmp)

read_array = np.array(read)             # stores raw data where rows mean indivisuals
N = len(read_array)
step_array = read_array.transpose()     # stores raw data where rows mean a given step/stage
S = len(step_array)

empirical = cvdiff(step_array)
plot_emp = []
for item in empirical:
    plot_emp.append([item])

deltas = []

for i in range(bootstraps):
    boot = bootrstrap(read_array)
    step_boot = boot.transpose()
    deltas.append(cvdiff(step_boot))

deltas_array = np.array(deltas)

if more_sig:
    t1, t2, t3 = significance(deltas_array)
else:
    t1 = significance(deltas_array)

results_directory = os.path.join(directory, 'Results')
if not os.path.exists(results_directory):
    os.mkdir(results_directory)

sns.barplot(plot_emp, color=('Grey'))
plt.boxplot(t1, positions = range(len(plot_emp)))

if more_sig:
    plt.boxplot(t2, positions = range(len(plot_emp)))
    plt.boxplot(t3, positions = range(len(plot_emp)))
    
plt.xlabel('Process step/stage')
plt.ylabel('CVabs - CVfrac')
plt.title(source_file)
plt.savefig(os.path.join(results_directory, source_file + '.png'))
plt.show()

sns.violinplot(deltas_array, color=('Orange'))
plt.xlabel('Process step/stage')
plt.ylabel('Simulated CVabs - CVfrac distributions')
plt.title(source_file + '_simul')
plt.savefig(os.path.join(results_directory, source_file + '_simul.png'))
plt.show()

# save bootstrap results to a csv file
if save_boots:
    with open(os.path.join(results_directory, source_file + '_bootsrapped.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for i in range(bootstraps):
            writer.writerow(deltas_array[i])

# save empirical significances results to a csv file
s = pval(empirical, deltas_array)
print(s)

if save_sigs:  
    with open(os.path.join(results_directory, source_file + '_significance.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(s)