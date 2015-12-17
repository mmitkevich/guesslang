#ifndef GUESSLANG_CHARSTAT_H__
#define GUESSLANG_CHARSTAT_H__

#include "deps.h"
#include <regex>
#include <limits>

namespace guesslang
{
    inline bool is_alpha(const QChar &c) {
        return QChar(c).isLetter();
    }

    template<typename V, typename T> V StaticCast(const T& d) {
        return static_cast<V>(d);
    }

    const double min_double = -1e50;

    template<typename T> bool True(const T& arg) { return true; }

//    template<typename V, typename T, typename L=QList<T>>
//    int max_index(const L& list,
//                  T* value=nullptr,
//                  std::function<V(const T&) transform=StaticCast<V,T>,
//                  const V& minval = std::numeric_limits<V>::min) {

//        int best = -1;
//        V best_val = minval;
//        for(int i=0; i<list.size(); i++) {
//            auto val = transform(list[i]);
//            if(best_val < val) {
//                best_val = val;
//                best = i;
//            }
//        }
//        return best;
//    }

    /// character probabilities
    class QCounter
    {
    public:
        typedef QHash<QString, int> Dict;
        typedef Dict::iterator iterator;
        typedef Dict::const_iterator const_iterator;

        /// counts character sequence of length m in input
        QCounter(int m = 1, std::function<bool(const QChar&)> filter = True<QChar>):
            total_(0),
            filter_(filter),
            m_(m)
        {
            assert(m>0);
        }

        /// probability of charcter
        double proba(QString c) {
            return double(cnt_.value(c, 0)) / total_;
        }

        int operator[](QString c) const {
            return cnt_[c];
        }

        QCounter& operator+=(const QCounter &rhs) {
            for(auto it = rhs.cnt_.cbegin(); it!=rhs.cnt_.cend(); it++) {
                cnt_[it.key()] += it.value();
            }
            return *this;
        }

        QCounter& operator-=(const QCounter &rhs) {
            for(auto it = rhs.cnt_.begin(); it!=rhs.cnt_.end(); it++) {
                cnt_[it.key()] -= it.value();
                if(cnt_[it.key()]==0) {
                    cnt_.remove(it.key());
                }
            }
            return *this;
        }

        void read_all(QTextStream& stream) {
            QChar c;
            while(!stream.atEnd()) {
                stream >> c;
                if(filter_(c))
                    append(c);
#ifndef NDEBUG
                else
                    qDebug() << "ignored " << c << endl;
#endif
            }
        }

        void append(QChar c) {
            buf_.push_back(c.unicode());
            if (buf_.size() >= m_) {
                QString str(m_, ' ');
                qCopy(buf_.begin(), buf_.end(), str.begin());
                buf_.pop_front();
                increment(str, 1);
            }
        }

        /// update counter
        void increment(const QString& s, int n) {
            //DEBUG("push " << s)
            total_ += n;
            cnt_[s] += n;
        }

        int total() const {
            return total_;
        }

        int size() const {
            return cnt_.size();
        }

        const_iterator cbegin() const {
            return cnt_.cbegin();
        }

        const_iterator cend() const {
            return cnt_.cend();
        }

        iterator begin() {
            return cnt_.begin();
        }

        iterator end() {
            return cnt_.end();
        }


    private:
        QQueue<QChar> buf_;
        int m_;
        Dict cnt_;
        int total_;
        std::function<bool(QChar)> filter_;
    };


    class QParagraph
    {
    public:
        QParagraph(int order=2, std::function<bool(const QChar&)> filter=True<QChar>, const QString& label=QString::null):
            order_(order),
            full_(order, filter),
            prefix_(order-1, filter),
            label_(label)
        {
            Q_ASSERT(order==2);
        }

        void append(QChar c) {
            full_.append(c);
            prefix_.append(c);
        }

        int read(QTextStream &stream, const QString &sep="\n", int min_chars=0, int max_chars=std::numeric_limits<int>::max()) {
            int n_chars = 0;
            while(n_chars<max_chars) {
                QChar c;
                if(stream.atEnd())
                    return n_chars;
                stream >> c;
                append(c);
                n_chars++;
                if(n_chars>=min_chars && sep.contains(c))
                    break;
            }
            return n_chars;
        }

        int alphabet_size() const {
            return prefix_.size();
        }

        const QString &label() const {
            return label_;
        }

        /// $Pr(w_1..w_ord|w_1..w_{ord-1})$
        double condi_logproba(QString f) {
            Q_ASSERT(f.size()==order_);
            auto m = prefix_.size();
            return (full_[f]+1)/(prefix_[f.left(1)]+m);
        }

        const QCounter& full() const {
            return full_;
        }

        const QCounter& prefix() const {
            return prefix_;
        }

        int order() {
            return order_;
        }

        QString str() const {
            return QString().sprintf("[f=%d,p=%d,m=%d]",full_.total(), prefix_.total(), prefix_.size());
        }

    protected:
        QCounter full_;
        QCounter prefix_;
        int order_;
        QString label_;
    };

    inline QTextStream & operator>>(QTextStream& stream, QParagraph& p)
    {
        for(QChar c: stream.readAll()) {
            p.append(c);
        }
        return stream;
    }

    const double min_likelihood = -1e50;

    /// language is set of paragraph samples
    ///
    class QLanguage : public QParagraph
    {
    public:
        typedef QList<QParagraph> List;
        typedef List::const_iterator const_iterator;

        QLanguage(int order = 2, const QString &label=QString()):
            QParagraph(order, True<QChar>, label),
            likelihood_(0)
        {

        }

        double likelihood() const {
            return likelihood_;
        }

        /// calculate probability that paragraph S belongs to language using order (j+1)
        /// $ln P(S|L) = \sum T(w_1..w_j,S)*\frac{T(w_1..w_j,L)+1}{T(w_1..w_{j-1},L)+m(L)}$
        /// m(L) = number of distinct characters in L
        /// j = j-gram length, e.g. for j=2 we analyze P(ab|a)
        double static calc_logproba(const QParagraph& s, const QList<QParagraph> &samples, int m) {
            double lprob = 0.;
            auto full = s.full();
            for(auto it=full.begin(); it!=full.end(); it++) {
                auto w = it.key();  // char sequence in sample of length j+1
                auto cnt = it.value();
                int full_cnt = 0, prefix_cnt = 0;
                for(const auto &ls: samples) {
                    full_cnt += ls.full()[w];
                    prefix_cnt += ls.prefix()[w.left(1)];
                }
                lprob += double(cnt) * double(full_cnt + 1) / double(prefix_cnt + m);
            }
            return lprob;
        }

        int static calc_alphabet_size(const QList<QParagraph>& samples1) {
            QCounter cnt;
            for(const auto &ls: samples1) {
                cnt += ls.prefix();
            }
            return cnt.size();
        }


        double static calc_likelihood(const QList<QParagraph>& samples1, int m) {
            double llh = 0.;//min_likelihood;
            //DEBUG("calc_likelihood(n_samples="<<samples1.size()<<")")
            if(samples1.size()>0) {
                //auto alphabet = union_alphabet(samples1);
                //auto m = alphabet.size();
                for(const auto &ls: samples1) {
                    llh += calc_logproba(ls, samples1, m);
                }
                /// here we should normalize: \sum e^{lp_k} = Z
                /// e^{lp1_k}=e^{lp_k}/Z
            }
            return llh;
        }


        /// append sample and recalculate likelihood
        /// \returns amount likelihood increased
        void append(const QParagraph &s, double *delta=nullptr) {
            samples_.append(s);
            prefix_ += s.prefix();
            full_ += s.full();
            if(delta) {
                likelihood_ += *delta;
            } else {
                likelihood_ = calc_likelihood(samples_, alphabet_size());
            }
        }

        /// remove sample and recalculate likelihood
        /// \param delta if specified, increment likelihood by that value. otherwise recalculate likelihood
        /// \returns sample
        QParagraph takeAt(int i, double *delta=nullptr) {
            QParagraph result = samples_.takeAt(i);
            prefix_-= result.prefix();
            full_ -= result.full();

            if(delta) {
                likelihood_ += *delta;
            } else {
                likelihood_ = calc_likelihood(samples_, alphabet_size());
            }
            return result;
        }

        double if_append(const QParagraph &s) const {
            QList<QParagraph> samples1(samples_);
            samples1.append(s);
            auto llh = calc_likelihood(samples1, calc_alphabet_size(samples1));
            return llh-likelihood_;
        }

        double if_removeAt(int i) const {
            QList<QParagraph> samples1(samples_);
            samples1.removeAt(i);
            auto llh = calc_likelihood(samples1, calc_alphabet_size(samples1));
            return llh-likelihood_;
        }

        int size() const {
            return samples_.size();
        }

        /// \returns index of what sample to remove to ge maximum likelihood increase
        int best_to_remove(double *delta = nullptr) const {
            int best = -1;
            double best_value = min_double;
            for(int i=0; i<samples_.size(); i++) {
                auto val = if_removeAt(i);
                if(val > best_value){
                    best = i;
                    best_value = val;
                }
            }
            if(delta)
                *delta = best_value;
            return best;
        }

        int best_to_remove2(double *delta = nullptr) const {
            int best = -1;
            double best_value = min_double;
            for(int i=0; i<samples_.size(); i++) {
                //auto val = if_removeAt(i);
                const auto val = -likelihood_for(samples_[i]);
                if(val > best_value){
                    best = i;
                    best_value = val;
                }
            }
            if(delta)
                *delta = best_value;
            return best;
        }

        double likelihood_for(const QParagraph &s, bool isnew=false) const{
            double lprob = 0.;
            QCounter f = full_;
            if(isnew)
                f+= s.full();
            QCounter p = prefix_;
            if(isnew)
                p+= s.prefix();
            auto m = p.size();
            for(auto it=s.full().cbegin(); it!=s.full().cend(); it++) {
                const auto &w = it.key();  // char sequence in sample of length j+1
                auto cnt = it.value();
                lprob += double(cnt) * double(f[w] + 1) / double(p[w] + m);
            }
            return lprob;
        }

        QList<QParagraph>::const_iterator cbegin() const {
            return samples_.cbegin();
        }

        QList<QParagraph>::const_iterator cend() const {
            return samples_.cend();
        }

        QString str() const {
            QString r;
            r.sprintf("L%s(ns=%3d,likely=%5g)\n", label().toUtf8().data(), size(), likelihood());
            QHash<QString, int> member_counts;
            for(const auto &s: samples_) {
                member_counts[s.label()] +=1;
            }
            for(decltype(member_counts)::iterator it = member_counts.begin(); it!= member_counts.end(); it++) {
                r.append(QString().sprintf(" %d %s|", it.value(), it.key().toUtf8().data()));
            }
            return r;
        }

    protected:
        List samples_;
        double likelihood_;
    };

    template<typename T> QTextStream& operator<<(QTextStream& stream, QLanguage &lang) {
        stream << lang.str();
        return stream;
    }


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
    class QClassifier
    {
    public:
        QClassifier(int order, int max_lang=50):
            likelihood_(0.),
            id_(0),
            n_samples_(0),
            max_lang_(max_lang),
            order_(order)

        {
            langs_.append(std::move(QLanguage(order, QString::number(id_++)))); // default language
        }

        double append(const QParagraph &s) {
            double delta = 0.;
            //int into_l = best_to_append(s, &delta);
            int into_l = 0;
            int min_size = 1000000;
            for(int l=0;l<langs_.size();l++){
                if(langs_[l].size()<min_size){
                    into_l=l;
                    min_size=langs_[l].size();
                }
            }
            langs_[into_l].append(s, &delta);
            likelihood_ += delta;
            if(langs_[into_l].size()==1 && langs_.size()<max_lang_) {   // always keep one empty unknown language with minimal likelihood for possible shuffling
                langs_.append(std::move(QLanguage(order_, QString::number(id_++))));
            }
            n_samples_++;
            //double metric = langs_[into_l].likelihood_for(s);
            //return metric;
            return delta;
        }

        int n_samples() const {
            return n_samples_;
        }

        int shuffle(int max_itr=1000, double rate=1e-3) {
            for(int itr=0; itr<max_itr; itr++) {
                int from_s = -1, into_s = -1;
                double delta = 0.;
                int from_l = best_to_remove(&delta, &from_s);
                DEBUG("rm from "<<langs_[from_l].str()<<"("<<langs_[from_l].size()<<") => "<<delta);
                QParagraph s = langs_[from_l].takeAt(from_s, &delta);
                int from_size = langs_[from_l].size();
                if(from_size==0) {
                    langs_.removeAt(from_l);        // remove empty
                    DEBUG("empty "<<from_l);
                }
                double change = delta;
                likelihood_ += delta;
                int into_l = best_to_append(s, &delta);
                DEBUG("ins into"<<langs_[into_l].str()<<" => "<<delta);
                langs_[into_l].append(s, &delta);
                likelihood_ += delta;
                change += delta;
                if(langs_[into_l].size()==1) {   // into was that empty class so append empty one
                    langs_.append(std::move(QLanguage(order_, QString::number(langs_.size()))));
                }
                double r = fabs(change/likelihood_);
                if(from_l!=into_l)
                    DEBUG("iter="<<itr<<" | LL="<<likelihood_<<" | n_lang="<<langs_.size()<<" | rate="<<r<<" | change="<<change<<" |from "<<from_l<<" |into "<<into_l);
                if(r<rate || change<0)
                    return itr;
            }
            return max_itr;
        }

        int size() {
            return langs_.size();
        }

        const QLanguage& operator[](int i) {
            return langs_[i];
        }

        double likelihood() {
            return likelihood_;
        }

        int best_to_remove(double *delta, int *best_sample) {
            int best = -1;

            double best_val = min_double;
            double d1 = 0.;
            for(int l=0; l<langs_.size(); l++) {
                d1 = min_double;
                int s = langs_[l].best_to_remove(&d1);
                if(d1>=best_val) {
                    best_val = d1;
                    best = l;
                    if(best_sample)
                        *best_sample = s;
                }
            }
            if(delta)
                *delta=best_val;
            return best;
        }

        /// \returns language to append s to get maximum likelihood increase
        int best_to_append(const QParagraph& s, double *delta) const {
//            return max_index<double>(langs_, delta, [s, langs_](const QLanguage& l){return l.if_append(s);});
            int best = -1;
            double best_val = min_double;
            for(int i=0; i<langs_.size(); i++) {
                //double val = langs_[i].if_append(s);
                double val = langs_[i].likelihood_for(s, true);
                //DEBUG("best_to_append "<<i<<"["<<langs_[i].size()<<"]"<<":"<<val<<":"<<langs_[i].likelihood());
                if(val>=best_val) {
                    best = i;
                    best_val = val;
                }
            }
            if(delta)
                *delta = best_val;
            //DEBUG("best_to_append->"<<best<<":"<<best_val);
            return best;
        }

        const QList<QLanguage> &languages() const {
            return langs_;
        }

        QString str() const {
            QString s;
            for(int i=0; i<langs_.size(); i++) {
                s += QString().sprintf("L%d:", i) + langs_[i].str()+"\n";
            }
            return s;
        }
    private:
        QList<QLanguage> langs_;
        double likelihood_;
        int id_;
        int n_samples_;
        int max_lang_;
        int order_;
    };
};

#endif //GUESSLANG_CHARSTAT_H__
