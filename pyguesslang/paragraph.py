#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import Counter
from itertools import chain
import re
import logging
import glob
import os

class Paragraph(object):
    def __init__(self, words, filename, seq=0):
        self.words = words
        self.filename = filename
        self.seq = seq

    def count(self, seqlen=1, alphabet=None):
        '''
        each word is treated as set of char sequences of seqlen characters, e.g.
            ^Hamlet$ == [^H, Ha, am, ml, le, t$]
        :param alphabet: all found characters should be in the alphabet specified (or ignored)
        :param seqlen: length of char sequence
        :return: dict(charseq->count)
        '''
        h = Counter()
        for w in self.words:
            #w = str([c for c in w0 if c in alphabet]) if alphabet else w0
            if len(w) > seqlen:
                for i in range(0, len(w)-seqlen):
                    syll = w[i:i+seqlen]
                    h[syll] += 1
        return h

    def __repr__(self):
        r = u"p %d file=%s words=%d -- %s..." % (self.seq, self.filename, len(self.words), ur' '.join(self.words[:10])[:60])
        return r.encode('utf-8')

    @staticmethod
    def loaddir(filemask="*.*", datadir=".", alphabet=set(), **kwargs):
        filenames = glob.glob(os.path.join(datadir, filemask))
        paragraphs = sorted(chain.from_iterable((Paragraph.loads(fn, alphabet, **kwargs) for fn in filenames)), key=lambda p:p.seq)
        return paragraphs

    @staticmethod
    def loads(filename, alphabet=set(), split_paragraph=ur"\n+", max_samples=None, min_chars=0, re_word=ur"\w+", exclude=ur"[\[\]{}\(\)\.\,\!]"):
        with open(filename) as fp:
            content = fp.read().decode('utf-8')
            #.lower()
            items = re.split(split_paragraph, content, re.UNICODE) if split_paragraph else [content]
            if max_samples:
                while len(items) > max_samples:
                    items = [items[2*i]+" "+items[2*i+1] for i in xrange(0, len(items)/2)]
            if min_chars>0:
                i = 0
                while i < len(items)-1:
                    if len(items[i]) < min_chars:
                        items[i] += u' '+items[i+1]
                        del items[i]
                    else:
                        i += 1

            paragraphs = []
            for text in items:
                def filter_bad_chars(w):
                    wr = ""
                    for c in w:
                        if c in exclude:
                            continue
                        wr += c
                    return wr
                words = [filter_bad_chars(w) for w in re.findall(re_word, text, re.UNICODE | re.DOTALL) if len(w) > 0]
                alphabet |= set(chain.from_iterable(words))
                text1 = u" ".join(words)
                paragraphs.append(Paragraph(words, filename, seq=len(paragraphs)))
            #logging.info("loaded %d paragraphs from %s", len(paragraphs), filename)
            return paragraphs



