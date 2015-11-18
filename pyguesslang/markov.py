from pyguesslang import *
import sklearn.cluster


def logprob(s, l, T, m):
    r = 0.0
    m1 = len(T[l][1])
    for w in T[s][0].iterkeys():
        p1 = float(T[l][0][w]+1)
        p2 = float(T[l][1][w[0:-1]]+m1)
        p_est = p1/p2
        lp_est = np.log(p_est)
        cnt = T[s][0][w]
        r += cnt * lp_est
    #print files[l], files[s], r
    return r

class MarkovChainClassifier(Classifier):
    def __init__(self, samples, alphabet, seqlen=2, to_label=None):
        self.alphabet = alphabet
        self.samples = samples
        self.N = len(samples)
        M = len(alphabet)
        self.seqlen = seqlen
        self.to_label = to_label or (lambda i, s: "S%d" % i)
        counts = [(pi.count(seqlen, self.alphabet), pi.count(seqlen-1, self.alphabet)) for pi in samples]
        self.dist_ls = np.zeros(shape=[self.N, self.N])
        for s, ps in enumerate(samples):
            for l, pl in enumerate(samples):
                lp = logprob(s, l, T=counts, m=M)
                self.dist_ls[s][l] = lp
            a = np.max(self.dist_ls[s])
            self.dist_ls[s] = -(self.dist_ls[s]-a)/a
        print self.dist_ls

    def fit(self, k, **kwargs):
        ql = Quality()
        self.dbs = sklearn.cluster.DBSCAN(eps=k)
        kls = self.dbs.fit_predict(self.dist_ls)
        for s in xrange(0, len(kls)):
            ql.inc(kls[s], self.to_label(s, self.samples[s]))
        # TODO: merge paragraphs in the samle klass and retrain

        return ql

    def test(self, samples):
        ql = Quality()
        #for i, (s, f) in enumerate(samples):
        #    sorted_dists = sorted(dists[i])
        #    ql["klass"] = kl = np.argmin(dists[i])
        #    ql["dist_rel_to_next"] = sorted_dists[0]/sorted_dists[1]
        #    ql.inc(kl, self.to_label(i, s))

        return ql

