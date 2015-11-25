import numpy as np
from scipy.special import gammaln
import math
import sys
import pylab as pl


def gamma_func_expansion_bb(x, n, a, b):
    # Beta binomial function for large values of x, n, a and b.
    # With large value of n and x, n choose x becomes "inf".
    # This is the gamma function expansion of the beta binomial distribution function.
    # from: stackoverflow.com/question/26935127/beta-binomial-function-in-python
    answer = gammaln(n + 1) + gammaln(x + a) + gammaln(n - x + b) + gammaln(a + b) - (
        gammaln(x + 1) + gammaln(n - x + 1) + gammaln(a) + gammaln(b) + gammaln(n + a + b))
    return math.exp(answer)


def extract_beta_binomial_parameters(group):
    # This function
    # * requires a "group" consisting of 16 samples, split into previous and latest measurement and then
    # divided into total and maligne isoforms.
    # * aggregates the measurements of day 1 to 4 to one sum (assuming the persons genes have not changed).
    # * removes empty and/or infinite entries.
    # * converts the input into the beta binomial parameters k, n, a and b.

    prev_m = group[:4]
    prev_t = group[8:12]
    latest_m = group[4:8]
    latest_t = group[12:16]

    # Remove samples of one day when m or t is inf
    prev_m = prev_m[prev_m != "inf"]
    prev_t = prev_t[prev_m != "inf"]

    prev_m = prev_m[prev_t != "inf"]
    prev_t = prev_t[prev_t != "inf"]

    latest_m = latest_m[latest_m != "inf"]
    latest_t = latest_t[latest_m != "inf"]

    latest_m = latest_m[latest_t != "inf"]
    latest_t = latest_t[latest_t != "inf"]

    # Remove empty cells and convert to float
    prev_m = prev_m[prev_m != ""].astype(float)
    prev_t = prev_t[prev_t != ""].astype(float)
    latest_m = latest_m[latest_m != ""].astype(float)
    latest_t = latest_t[latest_t != ""].astype(float)

    b = np.zeros((4))
    a = np.zeros((4))
    k = np.zeros((4))
    n = np.zeros((4))

    for i in range(4):
        b[i] = int(prev_t[i] - prev_m[i])
        n[i] = int(np.sum(latest_t[:]))
        k[i] = int(np.sum(latest_m[:]))
        a[i] = int(prev_m[i])

    # Aggregate all samples
    #prev_m = np.sum(prev_m)
    #prev_t = np.sum(prev_t)
    #latest_m = np.sum(latest_m)
    #latest_t = np.sum(latest_t)

    # Calculate Beta Binom Parameters and convert to int
    #b = int(prev_t - prev_m)
    #n = int(latest_t + latest_m)
    #k = int(latest_m)
    #a = int(prev_m)

    return k, n, a, b


def cdf_beta_binom(k, n, a, b):
    if k > 0 and n > 0 and a > 0 and b > 0:
        bb_mean = n * a / (a + b)

        cdf = 0
        bb_function = []
        for i in range(int(k)):
            cdf += gamma_func_expansion_bb(i, n, a, b)
            bb_function.append((i, gamma_func_expansion_bb(i, n, a, b)))
        bb_function = np.array(bb_function)
        cdf = min(cdf, 1.)

        pl.xlabel("Parameter k")
        pl.ylabel("beta binomial distribution")
        pl.plot(bb_function[:,0], bb_function[:,1])

        if bb_mean < int(k):
            probability = 1 - cdf
        else:
            probability = cdf

        bb_variance = n * a * b * (a + b + n) / math.pow((a + b), 2) / (a + b + 1)
        bb_stdev = math.sqrt(bb_variance)
        stdevs = abs(bb_mean - k) / bb_stdev
        print "Standard Deviations: %4.1f \tProbability %2.2e \tParameters: x=%i\tn=%i\ta=%i\tb=%i" % (
            stdevs, probability, k, n, a, b)
    else:
        print "not enough samples!"
        probability = 0
        stdevs = -1
    return probability, stdevs


def main():
    filename = "C:/Users/Hendrik/Downloads/exercise-master/exercise-master/raw_data.tsv"
    data = np.loadtxt(filename, dtype="S", delimiter="\t")
    header = data[0]
    data = data[1:, :]

    analyzed_data = data
    #analyzed_data = data[::10]
    probabilities = []
    standard_deviations = []

    for row in analyzed_data:
        k, n, a, b = extract_beta_binomial_parameters(row)
        for i in range(4):
            probability, stdevs = cdf_beta_binom(k[i], n[i], a[i], b[i])
        pl.legend(range(4))
        pl.show()
        probabilities.append(probability)
        standard_deviations.append(stdevs)

    probabilities = np.array(probabilities, dtype=float)
    standard_deviations = np.array(standard_deviations, dtype=float)
    person_id = np.arange(len(probabilities))

    results = np.vstack([person_id, probabilities, standard_deviations])
    results = np.transpose(results)
    results = results[results[:, 1].argsort()[::-1]]
    print results

    save_file = "C:/Users/Hendrik/Downloads/exercise-master/exercise-master/results.csv"
    np.savetxt(save_file, results, delimiter= ";")

    x = results[:, 1][results[:, 1] > 0]
    pl.plot(range(len(x)), x)
    pl.yscale("log")

    pl.figure()
    x = standard_deviations[standard_deviations > 0]
    pl.hist(x, bins=math.sqrt(len(x)))
    pl.show()


if __name__ == "__main__":
    sys.exit(main())