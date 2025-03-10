import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from time import time
from tqdm import tqdm
from scipy import stats
import math

# assumption: sampling with replacement!

def initialize_strands(n):
    """Initialize the number of copies for each strand to be constant- uniform distribution."""
    strand_counts = np.ones(n)
    return strand_counts

def create_graph(new_counts):
    """
    gets new_counts- new_counts[i] is the amount of time strand i was picked
    return x,y where:
    x is number of copies
    y is number of strands sampled x times
    """
    x = np.unique(new_counts)
    x = np.sort(x)

    # st = time()
    # y = np.zeros(x.shape[0])
    # for i in range(y.shape[0]):
    #     y[i] = np.sum(new_counts == x[i])
    # print(f"time regular {time()-st}")
    # st = time()

    unique_vals, counts = np.unique(new_counts, return_counts=True)
    y = np.zeros_like(x, dtype=int)  # Initialize result array

    # Match counts to corresponding values in x
    indices = np.searchsorted(unique_vals, x)
    mask = unique_vals[indices] == x  # Ensure only exact matches are counted
    y[mask] = counts[indices[mask]]

    # print(f"time otimized {time()-st}")
    # print(f"create_graph: {np.all(y_try==y)}")
    return x,y
    

def sample_strands(strand_counts, probabilities, M):
    """Take a sample of size M from the strands, based on their relative abundance."""

    sampled_strands = np.random.choice(len(strand_counts), size=M, p=probabilities)
    new_counts = np.bincount(sampled_strands, minlength=len(strand_counts))

    summed = strand_counts + new_counts

    return summed

def plot_distribution(x,y, iteration):
    """Plot the the requested graph"""
    plt.figure(figsize=(16, 10))

    plt.bar(x, y, color='blue', edgecolor='blue', alpha=0.7)
    plt.xlabel('Number of Copies')
    plt.ylabel('Number of Samples')
    plt.title(f'Samples per number of Copies - Iteration {iteration}')
    

    # Fit normal distribution to the data
    mu, std = stats.norm.fit(np.repeat(x, y))
    print(f"for iteration #{iteration}: mu = {mu} and std = {std}")
    
    # Generate points for the fitted normal distribution
    xmin, xmax = plt.xlim()
    x_fit = np.linspace(xmin, xmax, 100)
    y_fit = stats.norm.pdf(x_fit, mu, std) * sum(y)
    
    # Plot the fitted normal distribution
    plt.plot(x_fit, y_fit, 'r-', linewidth=2, label='Fitted Normal Distribution')


    plt.legend()
    plt.show()


def plot_multiple_distributions(mus, stds):
    plt.figure(figsize=(16, 10))
    x = np.linspace(min(mus) - 3*max(stds), max(mus) + 3*max(stds), 1000)
    
    for i, (mu, std) in enumerate(zip(mus, stds)):
        y = norm.pdf(x, mu, std)
        plt.plot(x, y, label=f'Distribution {i}')
        
        # Add index number next to the peak of each curve
        peak_x = mu
        peak_y = norm.pdf(peak_x, mu, std)
        plt.text(peak_x, peak_y, str(i), fontsize=12, verticalalignment='bottom')
    
    # plt.xlabel('X-axis')
    # plt.ylabel('Probability Density')
    # plt.title('Multiple Normal Distributions')
    plt.legend()
    plt.grid(True)
    plt.show()


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

    # plot_multiple_distributions(mus, stds)
    expected = calc_std_analytical(M, n, iterations)
    
    print(f"{expected=}")
    print(f"{stds=}")

    plot_arrays_for_index(stds, expected)

    # x,y = create_graph(strand_counts)
    # plot_distribution(x,y, iterations)

    
    # for i in tqdm(range(200)):
    #     probabs = strand_counts / ((i+1)*M)
            
    #     strand_counts = sample_strands(strand_counts, probabs, M)
    #     if i == 30:
    #         x,y = create_graph(strand_counts)
    #         plot_distribution(x,y, 30)

    
    # x,y = create_graph(strand_counts)
    # plot_distribution(x,y, 200)
    
    
if __name__ == "__main__":
    # Parameters
    n = int(1e3)      # Number of unique strands
    M = int(1e6)        # Sample size
    # print(M)
    iterations = 30    # Number of iterations

    simulate_sampling(n, M, iterations)
    