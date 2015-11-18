from pyguesslang import *
import numpy as np

if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    klassify(datadir="../data/accurat-corpus",testdir="../data/hamlet", seqlen=1, min_k=2, max_k=10, min_chars=100, max_samples=8)