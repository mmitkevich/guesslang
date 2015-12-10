#ifndef QCLASSIFIER_H
#define QCLASSIFIER_H

#include <functional>
#include "QQ/QQ.h"
#include <cstdlib>

const double min_double = -1e10;

const double min_likelihood = -1e50;

/// calculate probability that paragraph S belongs to language using order (j+1)
/// $ln P(S|L) = \sum T(w_1..w_j,S)*\frac{T(w_1..w_j,L)+1}{T(w_1..w_{j-1},L)+m(L)}$
/// m(L) = number of distinct characters in L
/// j = j-gram length, e.g. for j=2 we analyze P(ab|a)
double inline laplace_likelihood_estimate(const QParagraph& s, const QParagraph& l) {
    s.check();
    l.check();
    double lprob = 0.;
    int alphabet_size = l[1].size();
    if(alphabet_size==0)
        return min_likelihood;
    for(auto it=s[2].begin(); it!=s[2].end(); it++) { //: qqitems(s[2])
        auto w = it.key();  // 2-ngram
        auto cnt = it.value();
        lprob += double(cnt) * double(l[2][w] + 1) / double(l[1][w] + alphabet_size);
    }
#ifdef DEBUG4
    qqStdOut() << "laplace_estimate (" << s << "," << l << "," << alphabet_size << ")->"<<lprob<<"\n";
#endif
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
    QString centroid_label(int i) const {
        return centroid_labels_.value(i, QString("C%1") % i);
    }

    QString sample_label(int i, bool anon=true) const {
        return sample_labels_.value(i, anon ? QString("S") : QString("S%1")  % i);
    }

    const QVector<T>& centroids() const {
        return centroids_;
    }

    T& centroid(int i) {
        return centroids_[i];
    }

    const QVector<T>& samples() const {
        return samples_;
    }

    T& sample(int i) {
        return samples_[i];
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

    int sample_centroid(int is) const {
        return sample_centroid_[is];
    }

    int centroid_size(int ic) const {
        return centroid_size_[ic];
    }

    int n_samples() const {
        return samples_.size();
    }

    typedef QMap<int, QCounter<QString>> LabelCounter;

    LabelCounter labels_counter() const {
        QMap<int, QCounter<QString>> dict;

        for(int is=0; is<sample_centroid_.size(); is++) {
            int ic = sample_centroid_[is];
            QCounter<QString>& cnt = dict[ic];
            cnt.increment(sample_label(is), 1);
            Q_ASSERT(dict[ic][sample_label(is)]>0);
            Q_ASSERT(dict[ic].total()>0);
        }
        return dict;
    }

    QString str() const {
        QString str;

        auto lcdict = labels_counter();

        for(int ic=0; ic<centroids_.size(); ic++) {
            str << centroid_label(ic) << ": ";
            auto &lc = lcdict[ic];
            for(auto it = lc.begin(); it!=lc.end(); it++) {
                QString s = QString("%1 %2") % it.key() % it.value();
                str << s << " | ";
            }
            str << "\n";
        }

        return str;
    }

protected:
    QVector<T> centroids_;
    QMap<int,QString> centroid_labels_;
    QVector<T> samples_;
    QMap<int,QString> sample_labels_;
    QVector<int> sample_centroid_;
    QVector<int> centroid_size_;
};


template<typename T>
class QKMedoidsClassifier : public QClassifier<T>, public QQMatrix<double, QQVectorMatrixData<double>>
{
public:
    typedef QQMatrix<double, QQVectorMatrixData<double>> Matrix;
    typedef std::function<double(const T&, const T&)> Estimate;
    typedef QClassifier<T> Klass;

    using Klass::samples_;
    using Klass::sample_centroid_;
    using Klass::centroid_size_;
    using Klass::centroids_;
    using Klass::sample_labels_;
    using Klass::centroid;
    using Klass::sample;
    QKMedoidsClassifier(int max_centroids, const T& p, Estimate estimate, int max_iterations=10):
        QQMatrix(0, max_centroids),
        max_iterations_(max_iterations),
        min_convergence_rate_(1e-3),
        estimate_(estimate)
    {
        centroid_size_.fill(0, max_centroids);
        p.check();
        for(int i=0;i<max_centroids;i++)
            centroids_.append(p);
    }

    int append(QVector<T>&& vec, const QString& label = QString()) {
        int cnt = 0;
        for(T it: vec) {
            append(std::move(it), label);
            ++cnt;
        }
        return cnt;
    }

    /// insert paragraph into language maximizing likelihood
    int append(T &&sample, const QString& label = QString(), int c = -1) {
        //if(c==-1)
        //    c = rand() % nc();        // add to random centroid

        int is = nr();
        samples_.append(sample);
        sample_centroid_.resize(is+1);
        if(label.length())
            sample_labels_[is] = label;
        resize(is+1);
        //row(nc());
        double fit = min_likelihood;
        assign_centroid(is, c, fit);
        shuffle();
        return is;
    }

    Matrix calc_estimates_if_add(int is, int jc) {
        Q_ASSERT(jc<nc());
        Matrix ests(1, nc());   // l(is,kc);
        T &s = sample(is);
        for(int kc=0; kc<nc(); kc++)
        if(kc!=jc) {
            T c = centroid(kc);
            ests.set(0, kc, estimate_(s, c));
        }
        T c = centroid(jc);
        c.check();
        Q_ASSERT(c.order()==2);
        c +=s;
        c.check();
        ests.set(0, jc, estimate_(s, c));
#if 0
        qqStdOut() << QString("calc_est_if_add(%1,%2)") % is % jc
                   << "->" << ests <<"|\n";
        qqStdOut().flush();
#endif
        return ests;
    }

    Matrix find_best_estimates(int is, int &ic, double& fit) {
        Matrix best(1, nc());
//#if 1
        int jc; // target centroid
        int ic0 = ic;
        ic = rand() % nc(); // if nothing more useful
        double fit0 = fit;
        if(ic<0)
            fit = min_likelihood;
        for(jc=0; jc<nc(); jc++)
        if(jc!=ic0)
        {
            Matrix ests = calc_estimates_if_add(is, jc);
            auto max_est = ests.max();
            auto fit1 = ests[jc]-max_est;
#ifdef DEBUG4
            qqStdOut() << QString("if(s=%1)from(c=%2)into(c=%3)then(f=%4)willbe(f=%5)") % is % ic0 % jc % fit0 % fit1 << " | "<<ests<<" | oldsizes "<<dumps(centroid_size_)<<"\n";
            qqStdOut().flush();
#endif
            if(fit1-fit0>fit-fit0) {

                fit = fit1;
                ic = jc;
                best = ests;
            }
        }
//#endif
#ifdef DEBUG4
        qqStdOut() << QString("find_best_est(%1,%2,%3)") % is % ic % max_likelihood;
#endif
        return best;
    }

    void calc_estimates(int ic) {
        for(int is=0;is<nr();is++) {
            at(is,ic) = estimate_(sample(is),centroid(ic));
        }
    }

    bool assign_centroid(int is, int &ic, double &fit){
        //qqStdOut() << "assign_centroid("<<is<<","<<ic<<"), centroids: " << dumps(centroid_size_);
        int ic0 = ic;
        //Matrix ests = (ic<0) ? find_best_estimates(is, ic, fit) : calc_estimates_if_add(is, ic);
        Matrix ests = find_best_estimates(is, ic, fit);
        Q_ASSERT(ic>=0);
        Q_ASSERT(ic<nc());
        T c = centroid(ic);
        c.check();
        c += sample(is);
        c.check();
        setrow(is, ests);
        centroid_size_[ic]++;
        calc_estimates(ic); // updates estimates that changed
#if DEBUG2
        qqStdOut() << "MATRIX:\n" << (Matrix&)*this;
#endif
        // put sample into the priority queue, ordered by l_{c(s)}(s) - max_{c} l(s, c)
        samples_queue_.insert(fit, is);
        sample_centroid_[is] = ic;
        //qqStdOut() << QString("assign(fit=%1,s=%2,c=%3)") % fit % is % ic << "\n";
        return ic0!=ic;
    }

    double likelihood(int is, int ic) {
        return 0;
    }

    double likelihood() {
        return 0.;
    }

    /// shuffle samples between centroids and return if the shuffle increased likelihood
    int shuffle() {
        int ntotalshuffles=0;
        for(int ii=0; ii<max_iterations_; ii++) {
            int nshuffled = 0;
            //            qqStdOut() << "iter "<<ii;
            // take sample from the priority queue
            for(auto iq = samples_queue_.begin(); iq!=samples_queue_.end();) {
                //int is = samples_queue_.take(samples_queue_.firstKey());
                int is = iq.value();
                double fit = iq.key();
#ifdef DEBUG4
                qqStdOut() << "move "<<is <<"("<<fit<<")\n";
#endif
            //for(int is=0;is<samples_.size();is++) { // go and try moving all the samples
                /// to be changed:
                /// sample_centroid_[is]
                int ic = sample_centroid_[is];

                // remove from centroid & update changed likelihoods
                //T c;
                T c = centroid(ic);
                c.check();
                //T c = centroid(ic);
                iq = samples_queue_.erase(iq);

                const T &s = sample(is);
                c -= s;
                centroid_size_[ic]-=1;
                calc_estimates(ic);
                int ic0 = ic;
                double fit0 = fit;
                assign_centroid(is, ic, fit);    // assigns best centroid
                double rate = fit/fit0;
                if(ic!=ic0 && rate>1.0+1e-3) {
                    nshuffled++;
#ifdef DEBUG4
                    qqStdOut() << QString("shuffle(s=%1,c=%2->%3)") % is % ic0 % ic << "\n";
                    qqStdOut().flush();
#endif
                }

            }
            ntotalshuffles+=nshuffled;
            if(nshuffled==0)
                break;
        }
        return ntotalshuffles;
    }
private:
    int max_iterations_;
    double min_convergence_rate_;
    Estimate estimate_;
    QVector<double> sample_max_likelihood_;
    QVector<double> sample_likelihood_;
    QMultiMap<double, int> samples_queue_;
};


#endif // QCLASSIFIER_H
