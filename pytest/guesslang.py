import logging

from itertools import chain, groupby
import os
import numpy as np
from sklearn.cluster import KMeans

from pyguesslang import *


def counter_tostr(counter, fmt=" %5s:%4d", max_count=50, per_line=10):
    i = 0
    str=""
    for k, v in counter.most_common(max_count):
        str += fmt % (k, v)
        i += 1
        if i % per_line == 0:
            str += u"\n"
    return str

class Features(list):
    @staticmethod
    def tostr(samples, per_line=10, key=None, max_count=50):
        str = u"\n"
        for sample, features in sorted(samples, key=lambda (s, f): key(s)) if key else samples:
            str += repr(sample).decode('utf-8')
            str += counter_tostr(samples, per_line, max_count)
            str += u"\n"
        return str

    @staticmethod
    def histogram(samples):
        g = Counter()
        for s, f in samples:
            g += f
        print "converting counter"
        h = Features.counter_tohistogram(g)
        print "converted counter"
        return h

    @staticmethod
    def counter_tohistogram(counter):
        total = sum(counter.values())
        return {k: float(v)/total if total>0 else 0. for k, v in counter.iteritems()}

class CounterFeatures(Features):
    def __init__(self, counter, samples):
        super(CounterFeatures,self).__init__(map(counter, samples))

    def __repr__(self):
        fs = self.features()
        fstr = u" ".join(fs)
        str = "features count=%d:\n    %s...\n" % (len(fs), fstr[:100])
        str += CounterFeatures.tostr(self)
        return str

    def features(self):
        '''
        :return: sorted list of feature names
        '''
        return reduce(set.union, (f.keys() for s, f in self), set())

    def histogram(self):
        return Features.histogram(self)

    def to_features_array(self, fs=None):
        fs = fs or sorted(self.features())
        def histo(v):
            return Features.counter_tohistogram(v)
        lst = [[histo(v).get(f, 0) for f in fs] for s, v in self]
        arr = np.array(lst)
        return arr


def test_par(datadir="../data", testdir=None, seqlen=1, min_clusters=2, max_clusters=3):
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

    alphabet = set()
    paragraphs = Paragraph.loaddir(filemask="*.sample", datadir=datadir, alphabet=alphabet, max_samples=20, min_chars=None, split_paragraph=ur"\n")
    filenames = set((p.filename for p in paragraphs))
    print "files: %d\n%s" % (len(filenames), " ".join(filenames))
    print "alphabet: %s" % ur"".join(sorted(alphabet))
    print "total paragraphs %d" % len(paragraphs)
    samples = CounterFeatures(lambda s: (s, s.count(seqlen=seqlen, alphabet=alphabet)), paragraphs)
    #print u"%s" % Features.tostr(samples, key=lambda p:  p.seq)

    obs = samples.to_features_array()
    features = sorted(samples.features())
    #print obs

    kopt = 1
    ads=[]
    for k in xrange(min_clusters, max_clusters+1):
        if len(obs)<=k:
            break
        km = KMeans(n_clusters=k, verbose=False)
        ds = km.fit_transform(obs)
        cluster = {}
        name = {}

        for i in xrange(0, len(ds)):
            s, f = samples[i]
            kl = np.argmin(ds[i])
            #print ds[i], s.filename, kl  #, repr(s)
            cluster.setdefault(kl, Counter())[s.filename] += 1
        ads.append(np.sum(ds*ds)/k)
        print "**** K:%2d,\t\tN:%3d\t\tds^2: %10g (%3.2f)" % (k, len(obs), ads[-1], ads[-1]/ads[-2] if len(ads)>1 else 1.)
        print "\n".join(["K%d:%s" % (n, counter_tostr(c)) for n,c in cluster.items()])
        if testdir:
            test_paragraphs = groupby(Paragraph.loaddir(filemask="*.sample", datadir=testdir, alphabet=set(), max_samples=1, min_chars=None, split_paragraph=ur"\n"),
                                      lambda p:p.filename)
            for fn, pars in test_paragraphs:
                test_samples = CounterFeatures(lambda s: (s, s.count(seqlen=seqlen, alphabet=alphabet)), pars)
                test_obs = test_samples.to_features_array(sorted(samples.features()))
                ds = km.transform(test_obs)
                sds = sorted(ds[0])
                kl = np.argmin(ds)
                kn = sds[0]/sds[1]
                print "test samples %40s \tklass K%d dist=%5f=%5.2f*avg%5.2f*next" % (fn, kl, ds[0][kl], ds[0][kl]/(np.sum(ds[0])/k), kn)#, ds



if __name__ == "__main__":
    test_par(datadir="../data/accurat-corpus",testdir="../data/hamlet", seqlen=1, min_clusters=2, max_clusters=10)