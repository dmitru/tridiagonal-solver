
import matplotlib.pyplot as plt
from optparse import OptionParser

def plot(ys, xs):
    plt.axis([min(xs), max(xs), min(ys) - 1, max(ys) + 1])
    plt.plot(xs, ys)
    plt.show()
                
