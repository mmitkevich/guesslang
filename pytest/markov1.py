from sklearn.cluster import KMeans
from pyguesslang import *
import numpy as np
from itertools import imap
from sklearn.cluster import DBSCAN
import math


def test_markov(datadir="../data", testdir=None, seqlen=2, min_clusters=2, max_clusters=3):
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

    alphabet = set()
    paragraphs = Paragraph.loaddir(filemask="*.sample", datadir=datadir, alphabet=alphabet, max_samples=16, min_chars=None, split_paragraph=ur"\n")
    files = [pi.filename for pi in paragraphs]

    print "files: %d\n%s" % (len(files), " ".join(files))
    print "alphabet: %s" % ur"".join(sorted(alphabet))
    print "total paragraphs %d" % len(paragraphs)

    k = seqlen
    M = len(alphabet)
    N = len(paragraphs)

    counts = [(pi.count(k, alphabet), pi.count(k-1, alphabet)) for pi in paragraphs]

    def logprob(s, l, T, m):
        r = 0.0
        m1 =  len(T[l][1])
        for w in T[s][0].iterkeys():
            p1 = float(T[l][0][w]+1)
            p2 = float(T[l][1][w[0:-1]]+m1)
            p_est = p1/p2
            lp_est = math.log(p_est)
            cnt = T[s][0][w]
            r += cnt * lp_est
        #print files[l], files[s], r
        return r

    dist_ls = np.zeros(shape=[N, N])
    for s, ps in enumerate(paragraphs):
        for l, pl in enumerate(paragraphs):
            lp = logprob(s, l, T=counts, m=M)
            dist_ls[s][l] = lp
        a = np.max(dist_ls[s])
        dist_ls[s] = -(dist_ls[s]-a)/a

    print dist_ls

    ql = Quality()

    nneps = 50
    for neps in xrange(0, nneps):
        eps = 1.0/nneps*(neps+1)
        print "**** eps=%f"%eps
        dbs = DBSCAN(eps=eps)
        kls = dbs.fit_predict(dist_ls)
        for s in xrange(0, len(kls)):
            ql.inc(kls[s], files[s])

    print ql
if __name__ == "__main__":
    test_markov(datadir="../data/accurat-corpus", testdir="../data/hamlet", seqlen=1, min_clusters=2, max_clusters=10)