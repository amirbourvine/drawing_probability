import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from time import time
from tqdm import tqdm
from scipy import stats
import math

from sampling_script import create_graph,sample_strands

# assumption: sampling with replacement!
"""
This script is analytical vs experiments.
Plots the analytically calculated std values with the stds observed in the experiment.
"""


def plot_arrays_for_index(arr, expected):
    """
    Plots two arrays together with different colors and a legend.
    
    Parameters:
    arr (list or numpy array): First array to plot.
    expected (list or numpy array): Second array to plot.
    """
    plt.figure(figsize=(10, 5))
    plt.plot(range(len(arr)), arr, label='arr', color='blue', marker='o')
    plt.plot(range(len(expected)), expected, label='expected', color='red', marker='x')
    plt.xlabel("Index")
    plt.ylabel("Value")
    plt.title("Comparison of arr and expected")
    plt.legend()
    plt.grid(True)
    plt.show()


def calc_std_analytical(M, n, iterations):
    """
    gets M,n,iterations
    returns a list of the expected (analytically) std values for the iterations.
    """
    # new = old + old/i^2+1 + 2old/i
    var0 = M/n*(1-1/n)
    curr = 1

    variances = [curr]

    for i in range(1,iterations):
        curr = curr + curr/(i**2)+1+2*curr/i
        variances.append(curr)
    
    return [math.sqrt(v*var0) for v in variances]




def simulate_sampling(n, M, iterations):
    """Simulate the iterative sampling process."""
    
    probabs = np.ones(n) / n

    strand_counts = np.zeros(n)

    strand_counts = sample_strands(strand_counts, probabs, M)

    # print(new_counts.sum())
    
    mus = []
    stds = []

    for i in range(iterations):
        x,y = create_graph(strand_counts)
        # plot_distribution(x,y, i)

        mu, std = stats.norm.fit(np.repeat(x, y))
        print(f"for iteration #{i}: mu = {mu} and std = {std}")

        mus.append(mu)
        stds.append(std)

        probabs = strand_counts / ((i+1)*M)

        if( ((i+1)*M) != strand_counts.sum()):
            print("Error")
            exit(1)
        
            
        strand_counts = sample_strands(strand_counts, probabs, M)

    expected = calc_std_analytical(M, n, iterations)
    
    print(f"{expected=}")
    print(f"{stds=}")

    plot_arrays_for_index(stds, expected)



    
    
if __name__ == "__main__":
    # Parameters
    n = int(1e3)      # Number of unique strands
    M = int(1e6)        # Sample size
    # print(M)
    iterations = 30    # Number of iterations

    simulate_sampling(n, M, iterations)
    