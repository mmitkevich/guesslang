#ifndef QCLASSIFIER_H
#define QCLASSIFIER_H

#include <functional>
#include "QQ/QQ.h"
#include <cstdlib>
#include <cmath>

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
    const auto & l2 = l[2];
    const auto & l1 = l[1];
    const auto & s2 = s[2];
    for(auto it=s2.begin(); it!=s2.end(); it++) { //: qqitems(s[2])
        const auto &w = it.key();  // 2-ngram
        const auto &cnt = it.value();
        //const auto & cnt = s2[w];
        //auto cnt = s[2][w];
        auto c2 = l2.value(w, 0);
        auto c1 = l1.value(w[0], 0);
        auto lp1 = double(cnt) * double(c2 + 1) / double(c1 + alphabet_size);
        if(std::isnan(lp1)) {
            throw std::bad_exception();
        }
        lprob += lp1;
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
        return centroid_samples_[ic].size();
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
    QVector<QSet<int>> centroid_samples_;
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
    using Klass::centroid_samples_;
    using Klass::centroids_;
    using Klass::sample_labels_;
    using Klass::centroid;
    using Klass::sample;
    using Klass::sample_centroid;
    QKMedoidsClassifier(int max_centroids, const T& p, Estimate estimate, int max_iterations=10):
        QQMatrix(0, max_centroids),
        max_iterations_(max_iterations),
        min_convergence_rate_(1e-3),
        estimate_(estimate)
    {
        centroid_samples_.resize(max_centroids);
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
    int append(T &&s, const QString& label = QString()) {
        int is = nr();
        int ic = rand() % nc();

        samples_.append(s);
        sample_centroid_.append(ic);
        centroid_samples_[ic].insert(is);
        centroid(ic) += sample(is);

        if(label.length())
            sample_labels_[is] = label;

        resize(is+1);
        shuffle();
        return is;
    }

/*    Matrix calc_estimates_if_add(int ts, int jc, int is) {
        Q_ASSERT(jc<nc());
        Matrix ests(1, nc());   // l(is,kc);
        for(int kc=0; kc<nc(); kc++)
        if(kc!=jc) {
            T c = centroid(kc);
            ests.set(0, kc, estimate_(sample(ts), c));
        }
        T c = centroid(jc);
        c.check();
        Q_ASSERT(c.order()==2);
        c += sample(is);
        c.check();
        ests.set(0, jc, estimate_(sample(ts), c));
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
            double fit1 = 0.;
            Matrix ests(1, nc());   // l(is,kc);
            for(int kc=0; kc<nc(); kc++)
            if(kc!=jc) {
                T c = centroid(kc);
                ests.set(0, kc, estimate_(sample(ts), c));
            }
            T c = centroid(jc);
            c.check();
            Q_ASSERT(c.order()==2);
            c += sample(is);
            c.check();

            ests.set(0, jc, estimate_(sample(ts), c));

            for(int ts=0; ts<nr(); ts++)
            if(js!=is && sample_centroid(js)==ic0) {
                Matrix ests = calc_estimates_if_add(ts, jc, is);
                auto max_est = ests.max();
                fit1 += ests[jc]-max_est;
            }
            // and all other samples affected by addition
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
*/
    double update_estimates(Matrix& ests, int is, int kc, double e) {
        auto e0 = normalized_estimates(*this, is, kc);
        ests.set(is, kc, e);
        auto e1 = normalized_estimates(ests, is, kc);
        return e1-e0;
    }

    double normalized_estimates(Matrix &ests, int is, int kc) {
        auto em = ests.row(is).max();
        auto ne0 = ests.at(is, kc);
        auto ne = ne0 - em;
        if(std::isnan(ne)) {
            qqStdOut()<<em<<ne;
        }

        return ne;
    }

    int takefrom_centroid(double &dfitr) {
        int is_best = -1;
        Matrix best = *this;
        dfitr = 0.;
        is_best = nr()-1;
        for(int is=0; is<nr(); is++) {
            int ic = sample_centroid(is);
            Q_ASSERT(ic!=-1);
            // removed from ic affects all in ic
            auto c = centroid(ic);
            auto& s = sample(is);
            auto cminus = c - s;//centroid(ic0)-sample(is);
            double dfit = 0.;
            Matrix ests = *this;        // let's calc new ests
            for(int js=0; js<nr(); js++) {  // update all samples in the cluster
                if(sample_centroid(js) == ic) {
                    auto e = estimate_(sample(js), cminus);
                    auto f = update_estimates(ests, js, ic, e);
                    dfit += f;
                }
            }
            if(dfit>dfitr) {
                dfitr = dfit;
                is_best = is;
            }
        }

        if(is_best!=-1) {
            copy_from(best);
            int ic_best = sample_centroid(is_best);
            sample_centroid_[is_best] = -1;
            centroid_samples_[ic_best].remove(is_best);
            centroid(ic_best) -= sample(is_best);

            qqStdOut() << QString("takefrom(S%1, L%2, dfit=%3)") % is_best % ic_best % dfitr << "\n";
        }
        return is_best;
    }

    /// return new centroid assigned
    /// updates dfiti to increment in the fitness
    int assign_centroid(int is, double &dfiti){
        Q_ASSERT(sample_centroid(is)==-1);
        Matrix best = *this;
        auto& s = sample(is);
        dfiti = min_likelihood;
        int ic_best = -1;

        for(int jc=0; jc<nc(); jc++)    // look through target centroids
        //        if(jc!=ic0)                 // different from current
        {
            Matrix ests = *this;        // let's calc new ests
            double dfit = 0.;

            auto cplus = centroid(jc)+sample(is);

            // appended into jc
            auto e = estimate_(sample(is), cplus);
            auto f = update_estimates(ests, is, jc, e);    // added to target

            dfit += f;

            // new sample affects estimates for all samples in the target centroid
            for(int js=0; js<nr(); js++) {
                if(sample_centroid(js)==jc) {
                    e = estimate_(sample(js), cplus);
                    f = update_estimates(ests, js, jc, e);
                    dfit += f;
                }
            }

            if(ic_best<0 || dfit>dfiti) {
                dfiti = dfit;
                ic_best = jc;
                best = ests;
            }
        }
        if(ic_best!=-1) {
            copy_from(best);
            sample_centroid_[is] = ic_best;
            centroid_samples_[ic_best].insert(is);
            centroid(ic_best) += s;
        }


#if DEBUG2
        qqStdOut() << "MATRIX:\n" << (Matrix&)*this;
#endif
        qqStdOut() << QString("assign(S%1, L%2, dfit=%3)") % is % ic_best % dfiti << "\n";
        static int flushcnt = 0;
        if(++flushcnt>1000) {
            qqStdOut().flush();
            flushcnt=0;
        }

        return ic_best;
    }



    double likelihood(int is, int ic) {
        return 0;
    }

    double likelihood() {
        return 0.;
    }

    double whatif(Matrix& ests, int ic, T& cwhatif, int is_skip=-1) {
        auto dfit = 0.;
        foreach(int js, centroid_samples_[ic]) {
            if(js!=is_skip) {
                Q_ASSERT(sample_centroid(js) == ic);
                auto e = estimate_(sample(js), cwhatif);
                auto f = update_estimates(ests, js, ic, e);
                if(std::isnan(f)) {
                    qqStdOut()<<f;
                }

                dfit += f;
            }
        }
        if(std::isnan(dfit)) {
            qqStdOut()<<dfit;
        }

        return dfit;
    }

    /// shuffle samples between centroids and return if the shuffle increased likelihood
    int shuffle() {
        int nshuff = 0;
        static int nshuff_total = 0;
        auto dfit_best = 0.;
        auto dfit_prev = 0.;
        int ic_best = -1;
        int is_best = -1;
        //qqStdOut()<<"shuffle("<<max_iterations_<<")\n";
        for(int ii=0; ii<max_iterations_; ii++) {
            Matrix best = *this;
            int jc = rand() % nc();
            //for(int jc=0; jc<nc(); jc++)    // look through target centroids
            //if(jc!=ic)                 // different from current
            {
                //for(int is=0; is<nr(); is++)
                int is = rand() % nr();
                {
                    auto& s = sample(is);
                    auto ic = sample_centroid(is);
                    auto cminus = centroid(ic) - s;//centroid(ic0)-sample(is);
                    Matrix ests = *this;        // let's calc new ests
                    auto dfitr = whatif(ests, ic, cminus, is);


                    auto cplus = centroid(jc)+s;
                    // appended into jc
                    auto e = estimate_(sample(is), cplus);
                    auto dfiti = update_estimates(ests, is, jc, e);    // added to target
                    if(std::isnan(dfiti)) {
                        qqStdOut()<<dfiti;
                    }

                    dfiti += whatif(ests, jc, cplus);


                    if(ic_best<0 || dfiti+dfitr>dfit_best) {
                        dfit_best = dfiti+dfitr;
                        if(std::isnan(dfit_best)) {
                            qqStdOut()<<dfiti<<","<<dfitr;
                        }
                        ic_best = jc;
                        is_best = is;
                        best = ests;
                    }
                } // for jc...
            } // for is0...
            if(ic_best!=-1 && ic_best!=sample_centroid(is_best)) {
                copy_from(best);
                auto ic_from = sample_centroid(is_best);
                centroid(ic_from) -= sample(is_best);
                centroid_samples_[ic_from].remove(is_best);

                sample_centroid_[is_best] = ic_best;
                centroid_samples_[ic_best].insert(is_best);
                centroid(ic_best) += sample(is_best);

                nshuff++;
                if(nshuff_total++ % 1==0) {
                    qqStdOut()<<QString("[%5] shuffle %1 from %2 into %3 profit %4\n") % is_best % ic_from % ic_best % dfit_best % nshuff_total;
                     qqStdOut().flush();
                }

            }

            auto ratio = dfit_prev==0?0.:dfit_best/dfit_prev-1.;
            if(dfit_prev>0 && dfit_best>0 && ratio>0 && ratio<1e-3) {
                qqStdOut()<<QString("[%1] %2")%nshuff_total%ratio;
                return nshuff;
            }
            dfit_prev = dfit_best;
        } // for ii...

        return nshuff;
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
