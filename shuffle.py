import numpy as np
import math
import sys
import pylab as pl

def main():
    filename = "C:/Users/Hendrik/Downloads/exercise-master/exercise-master/raw_data.tsv"
    data = np.loadtxt(filename, dtype="S", delimiter="\t")
    header = data[0]
    data = data[1:, :]

    analyzed_data = data
    #analyzed_data = data[::10]
    findings = []

    for row in analyzed_data:
        if "" in row:
            print "missing value"
        elif "inf" in row:
            print "infinite value"
        else:
            row = row.astype(float)
            august_ratio = row[:4]/row[8:12]
            august_mean = np.mean(august_ratio)
            december_ratio = row[4:8]/row[12:16]
            december_mean = np.mean(december_ratio)
            both_ratios = np.hstack([august_ratio, december_ratio])

            samples = 10000
            sample_mean_diff = np.zeros(samples)
            for i in range(samples):
                random_sample = np.zeros((4))
                for j in range(4):
                    random_sample[j] = both_ratios[np.random.randint(8)]
                sample_mean_diff[i] = np.mean(random_sample-august_mean)
            #print august_mean, december_mean, december_mean-august_mean
            #print sample_mean_diff
            larger_then = np.sum(sample_mean_diff >= december_mean-august_mean)*1./samples
            smaller_then = np.sum(sample_mean_diff <= december_mean-august_mean)*1./samples
            #print larger_then, smaller_then

            significance = min(larger_then, smaller_then)
            print significance
            findings.append(significance)
            #if significance >= .05:
                #print larger_then, smaller_then, sample_mean_diff
                #pl.hist(sample_mean_diff)
                #pl.show()

    findings = np.array(findings, dtype=float)
    print findings, findings.argsort()
    print "Mean significance %2.2f" % np.mean(findings)
    print "Significant samples %i" % np.sum(findings>0.05)
    print "Fraction of significant samples %2.2f" % np.mean(findings>0.05)
    pl.plot(range(len(findings)), findings[findings.argsort()[::-1]])
    pl.figure()
    pl.hist(findings)
    pl.show()







    return

if __name__ == "__main__":
    sys.exit(main())