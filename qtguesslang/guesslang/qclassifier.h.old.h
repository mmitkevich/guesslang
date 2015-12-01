#ifndef QCLASSIFIER_H
#define QCLASSIFIER_H


#include "QQ/qqmatrix.h"

const double min_double = -1e10;

const double min_likelihood = -1e50;

/// calculate probability that paragraph S belongs to language using order (j+1)
/// $ln P(S|L) = \sum T(w_1..w_j,S)*\frac{T(w_1..w_j,L)+1}{T(w_1..w_{j-1},L)+m(L)}$
/// m(L) = number of distinct characters in L
/// j = j-gram length, e.g. for j=2 we analyze P(ab|a)
template<typename P, typename L>
double inline laplace_likelihood_estimate(P& s, L& l, int alphabet_size) {
    double lprob = 0.;
    if(alphabet_size==0)
        return min_likelihood;
    for(auto it: items(s[2])) {
        auto w = it.key();  // 2-ngram
        auto cnt = it.value();
        lprob += double(cnt) * double(l[2][w] + 1) / double(l[1][w] + alphabet_size);
    }
    return lprob;
}


typedef bool (*QCharPredicate)(const QChar&);


/// step t. we have K centroids (languages) $\{L_k\}, k=1..K$
/// each centroid $L_k$ contains samples $S^k_i, i=1..I_k$
/// step t+1. we try to maximize likelihood of clustering. So we move samples around languages.
/// For each sample S in each language L we have probability P(S|L).
/// let S_1 = argmin_{S,L} P(S|L) be sample with minimal probability.
/// move S_1 into language L_1 = argmax_{L} P(S1|L_1)
/// now we have new likelihood $ln PL = \sum_{L} \sum_{S \in L} ln P(S|L)$
/// we repeat until that likelihood stops to converge to its maximum
/// There could be situation that we move last sample from some language k to another language. That would mean we decrease K->K-1
/// We should have possibility to move each sample
template<typename T>
class QClassifier
{
public:
    QClassifier()
    {
    }

    QString centroid_label(int i) const {
        return centroid_labels_.value(i, QString("L%1") % i);
    }

    QString sample_label(int i) const {
        return sample_labels_.value(i, QString("S%1")  % i);
    }

    const QVector<T>& centroids() const {
        return centroids_;
    }

    const QVector<T>& samples() const {
        return samples_;
    }

    /// distance from sample to centroid
    ///virtual double operator()(int s, int c) = 0;

//        double likelihood() {
//            double val =  reduce<double>(langs_, [&](const QString &lbl, const QLanguage &lang, double &result)  {
//                result += lang.likelihood();
//                return true;
//            });
//            return val;
//        }

//        double likelihood_if(QSet<QLanguage&> removed, QSet<QLanguage&> inserted) {
//            double val =  reduce<double>(langs_, [&](const QString &lbl, const QLanguage &lang, double &result)  {
//                if(!removed.contains(lang))
//                    result += lang.likelihood();
//                else
//                    result += replaced.likelihood();
//                return true;
//            });
//            return val;
//        }

//        double likelihood_filtered(std::function<bool(const QLanguage&)> filter) {
//            double val =  reduce<double>(langs_, [&](const QString &lbl, const QLanguage &lang, double &result)  {
//                if(filter(lang.label))
//                    result += lang.likelihood();
//                return true;
//            });
//            return val
//        }


    QString str() const {
        return dumps(centroids_);
    }

protected:
    QVector<QParagraph> centroids_;
    QMap<int,QString> centroid_labels_;
    QVector<QParagraph> samples_;
    QMap<int,QString> sample_labels_;
};


template<typename T>
class QKMedoidsClassifier : public QClassifier<T>//, public QQMatrix<double>
{
public:
    QKMedoidsClassifier(int n_centroids) : QQMatrix(0, n_centroids) {
        max_iterations_ = 1000;
        min_convergence_rate_ = 1e-3;
    }

    template<typename D> int append(D& vec) {
        int cnt = 0;
        for(auto it: vec) {
            append(it);
            ++cnt;
        }
        return cnt;
    }

    /// insert paragraph into language maximizing likelihood
    int append(const T &&sample) {
//            double best = -1e10;
//            QLanguage lang_best;
//            double factor = -1e10;
//            for(auto it=langs_.begin();it!=langs_.end();it++) {
//                QLanguage lang1 = it.value() + s;

//                // log s
//                // log(s+e^l) = log s+ log(1+e^l/s)
//                auto est = lang1.likelihood(s);
//                factor += log(1.+exp(est-factor));
//                auto l_if_repl = likelihood_if_replace(lang1);
//                auto l_now = likelihood();
//                auto delta = l_if_repl - l_now;
//                if(delta>best) {
//                    lang_best = lang1;
//                }
//            }
//            langs_[lang_best.label()] = lang_best;
//            DEBUG(" replaced "<<langs_[lang_best.label()].str());
//            n_samples_++;
//            //DEBUG("ns "<<n_samples_<<" "<<s.str());
        //return *this;
        return 0;
    }

    /// return likelihood (log-probability) that sample s belongs to cluster c
    double likelihood(int s, int c) {
        return (*this)(s, c);
    }

    double likelihood() {
        return 0.;
    }

    int shuffle() {
//            for(int itr=0; itr<max_itr; itr++) {
//                double best = -1e10;
//                QLanguage from_best, into_best;
//                QParagraph sample_best;
//                double start = likelihood();
//                for(auto ifrom=langs_.begin();ifrom!=langs_.end();ifrom++)
//                    for(auto is: ifrom.value().samples()){
//                        double factor = 0.;
//                        for(auto iinto=langs_.begin();iinto!=langs_.end();iinto++) {
//                            QLanguage into1 = iinto.value() + is;
//                            QLanguage from1 = ifrom.value() - is;

//                            auto l_if_repl = likelihood_if_replace(into1, from1);
//                            auto l_now = likelihood();
//                            auto delta = l_if_repl - l_now;
//                            if(delta>best) {
//                                from_best = from1;
//                                into_best = into1;
//                                sample_best = is;
//                                best = delta;
//                            }
//                        }
//                   }
//                langs_[from_best.label()] = from_best;
//                langs_[into_best.label()] = into_best;
//                calc_factor(sample_best);
//                double rate = best/likelihood();
//                DEBUG(" moved "<<sample_best.str()<<" from "<<from_best.str()<<" into "<<into_best.str()<<" rate "<<QString::number(rate));
//                if(fabs(rate)<minrate)
//                    break;
//            }
        return 0;
    }
private:
    int max_iterations_;
    double min_convergence_rate_;
};


#endif // QCLASSIFIER_H
