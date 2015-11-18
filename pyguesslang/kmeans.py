from pyguesslang import *
import sklearn.cluster

class KMeansClassifier(Classifier):
    def __init__(self, paragraphs, alphabet, seqlen=1, to_label=None, features=None):
        samples = Samples(lambda s: (s, s.count(seqlen=seqlen, alphabet=alphabet)), paragraphs)
        self.alphabet = alphabet
        self.samples = samples
        self.N = len(samples)
        self.M = len(alphabet)
        self.to_label = to_label or (lambda i, s: "S%d" % i)
        self.features = features or sorted(samples.features())
        self.features_array = samples.to_features_array()

    def fit(self, k, **kwargs):
        ql = Quality()
        self.km = sklearn.cluster.KMeans(n_clusters=k, verbose=False)
        dists = self.km.fit_transform(self.features_array)
        for i, (s, f) in enumerate(self.samples):
            kl = np.argmin(dists[i])
            #print ds[i], s.filename, kl  #, repr(s)
            ql.inc(kl, self.to_label(i, s))
        ql.fitness = np.sum(dists*dists)/k
        return ql

    def test(self, samples):
        ql = Quality()
        farr = samples.to_features_array(self.features)
        dists = self.km.transform(farr)
        for i, (s, f) in enumerate(samples):
            sorted_dists = sorted(dists[i])
            ql["klass"] = kl = np.argmin(dists[i])
            ql["dist_rel_to_next"] = sorted_dists[0]/sorted_dists[1]
            ql.inc(kl, self.to_label(i, s))
        return ql
