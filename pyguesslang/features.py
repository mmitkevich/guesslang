import numpy as np
from itertools import chain, groupby
from collections import Counter
from pyguesslang import *

def features_tostr(counter, fmt=" %5s:%4d (%4.2f)", max_count=50, per_line=10):
    i = 0
    sm = sum(counter.values())
    str = "%4d --" % (sm)
    for k, v in (counter.most_common(max_count) if hasattr(counter,'most_common') else counter.iteritems()):
        str += fmt % (k, v, 1.0*v/sm)
        i += 1
        if i % per_line == 0:
            str += u"\n"
        if i>max_count:
            break
    return str

class Samples(list):
    @staticmethod
    def tostr(samples, per_line=10, key=None, max_count=50):
        str = u"\n"
        for sample, features in sorted(samples, key=lambda (s, f): key(s)) if key else samples:
            str += repr(sample).decode('utf-8')
            str += features_tostr(samples, per_line, max_count)
            str += u"\n"
        return str

    @staticmethod
    def histogram(samples):
        g = Counter()
        for s, f in samples:
            g += f
        print "converting counter"
        h = Samples.tohisto(g)
        print "converted counter"
        return h

    @staticmethod
    def tohisto(counter):
        total = sum(counter.values())
        return {k: float(v)/total if total>0 else 0. for k, v in counter.iteritems()}

    def __init__(self, extract_features, samples):
        super(Samples, self).__init__(map(extract_features, samples))

    def __repr__(self):
        fs = self.features()
        fstr = u" ".join(fs)
        str = "features count=%d:\n    %s...\n" % (len(fs), fstr[:100])
        str += Samples.tostr(self)
        return str

    def features(self):
        '''
        :return: sorted list of feature names
        '''
        return reduce(set.union, (f.keys() for s, f in self), set())

    def histogram(self):
        return Samples.histogram(self)

    def to_features_array(self, fs=None, norm=None):
        fs = fs or sorted(self.features())
        norm = norm or Samples.tohisto
        lst = [[norm(v).get(f, 0) for f in fs] for s, v in self]
        arr = np.array(lst)
        return arr


class Quality(dict):
    def __init__(self):
        self.klasses = {}
        self.fitness = -1e100

    def inc(self, k, label):
        self.klasses.setdefault(k, Counter())[label] += 1
        #print label, "K%d" % k

    def label(self, k):
        return self.klasses[k].most_common(1)[0][0]

    def klass(self, label):
        for k in self.klasses.keys():
            if self.label(k) == label:
                return k
        return None

    def __repr__(self):
        s1 = "\n".join(["K%d -- %s" % (k, features_tostr(c)) for k, c in self.klasses.iteritems()])
        #s2 = ",".join(["%s=%r" % (k,v) for k,v in self.iteritems()])
        #if len(s2) > 0:
        #    s1 +="\n("+s2+")"
        return s1

    def labels(self):
        return reduce(set.union, (set(c.keys()) for c in self.klasses.values()),set())

    def cleanness(self):
        fls = 0.0
        for kl, c in self.klasses.iteritems():
           return
    def report(self, vql):
        s = ""
        for lbl in self.labels():
            kl = self.klass(lbl)
            c = vql.klasses.get(kl)
            if c:
                s += "%s -- %s\n" % (lbl, features_tostr(c) if c else "n/a")
        return s

class Classifier(object):
    def __init__(self, data, alphabet, to_label=None):
        raise NotImplementedError()

    def fit(self, k, **kwargs):
        raise NotImplementedError()

    def test(self, samples):
        raise NotImplementedError()
