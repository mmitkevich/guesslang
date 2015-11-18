from pyguesslang import *
import logging

def klassify(klassifier, datadir="../data", testdir=None, seqlen=1, min_k=2, max_k=10, step_k=1, **kwargs):
    alphabet = set()
    paragraphs = Paragraph.loaddir(filemask="*.sample", datadir=datadir, alphabet=alphabet, **kwargs)
    filenames = set((p.filename for p in paragraphs))

    N = len(paragraphs)
    F = len(filenames)

    def strip_fn(s):
        items = re.split(r"[-\.]", os.path.basename(s))
        return items[0]+"-"+items[1] if len(items) > 1 else items[0]

    def to_label(i, s):
        return strip_fn(s.filename)

    print "input files: %d\n%s" % (F, " ".join(map(strip_fn, filenames)))
    print "alphabet: %s" % ur"".join(sorted(alphabet))
    print "total paragraphs %d" % N

    klsf = klassifier(paragraphs, alphabet, seqlen=seqlen, to_label=to_label)

    fits = []

    k = min_k
    while k < max_k:
        ql = klsf.fit(k=k)
        fits.append(ql.fitness)
        print "**** TRAIN: K = %3g, N = %3d, F = %3d, fit[i] = %10g = %3.2f*fit[i-1]" % (k, N, F, fits[-1], fits[-1]/fits[-2] if len(fits)>1 else 1.)
        print repr(ql)

        if testdir:
            test_paragraphs = Paragraph.loaddir(filemask="*.sample", datadir=testdir, alphabet=alphabet, **kwargs)
            V = len(test_paragraphs)
            test_samples = Samples(lambda s: (s, s.count(seqlen=seqlen, alphabet=alphabet)), test_paragraphs)
            tql = klsf.test(test_samples)
            print "**** CROSS: K=%3g, V=%3d" % (k, V)
            print ql.report(tql)

        k += step_k
