from sklearn.cluster import KMeans
from pyguesslang import *
import numpy as np
from itertools import imap
from sklearn.cluster import DBSCAN
import math


if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    klassify(MarkovChainClassifier, datadir="../data/accurat-corpus",testdir="../data/hamlet", seqlen=1, min_k=0.01, max_k=0.7, step_k=0.01, min_chars=300, max_samples=50)