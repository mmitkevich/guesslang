#ifndef QCOUNTER_H
#define QCOUNTER_H

#include <QtCore/QMap>

/// character probabilities
template<typename T>
class QCounter : public QHash<T, int>
{
public:
    //typedef QHash<QString, int> Dict;
    //typedef Dict::iterator iterator;
    //typedef Dict::const_iterator const_iterator;
    typedef QHash<T,int> Map;

    /// counts character sequences
    QCounter(int order=-1):
        total_(0),
        order_(order)
    {
    }

    /*QCounter(const QCounter &other):
        Map(other),
        order_(other.order_)
    {
        Q_ASSERT(order_==other.order_);
        other.check();
        check();
    }

    QCounter(QCounter &&other) {
        Q_ASSERT(order_==other.order_);
        other.check();
        swap(other);
    }*/

    void swap(QCounter &other) {
        Q_ASSERT(order_==other.order_);
        other.check();
        Map::swap(other);
        qSwap(total_, other.total_);
        qSwap(order_, other.order_);
        check();
    }

    void check() const{
#ifdef DEBUG3
        for(auto it=Map::begin();it!=Map::end();it++){
            T k = it.key();
            if(order_!=k.size()) {
                Q_ASSERT(0);
            }
        }
#endif
        if(total_<0)
            throw std::bad_exception();

    }

/*    QCounter &operator=(const QCounter &other) {
        Q_ASSERT(order_==other.order_);
        other.check();
        check();
        Map::operator =(other);
        order_ = other.order_;
        check();
        return *this;
    }
*/
    /// probability of charcter sequence
    double probability(const T &s) {
        return double(value(s, 0)) / total_;
    }

    /// add other counter into this counter
    QCounter& operator+=(const QCounter &other) {
        check();
        other.check();
        for(auto it=other.begin(); it!=other.end(); it++) {
            const auto &k = it.key();
            const auto &v = it.value();
            increment(k, v);
        }
        check();
        return *this;
    }


    QCounter operator+(const QCounter &other) const {
        check();
        other.check();
        QCounter result = *this;
        result.check();
        result+=other;
        result.check();
        return result;
    }

    /// subtract rhs counter from this counter
    QCounter& operator-=(const QCounter &other) {
        other.check();
        check();
        for(auto it=other.begin(); it!=other.end(); it++) {
            auto k = it.key();
            auto v = it.value();
            increment(k, -v);
        }
        check();
        return *this;
    }

    QCounter operator-(const QCounter &other) const {
        QCounter result = *this;
        result-=other;
        return result;
    }

    /// increment counter for sequence specified
    int increment(const T& w, const int n=1) {
        //DEBUG("inc " << dumps(*this));
        total_ += n;
        int a = Map::value(w, 0);
#ifdef DEBUG3
        check();
        if(order_!=-1) {
            Q_ASSERT(order_==w.length());
            if(Map::size()>0) {
                const T &kk = Map::firstKey();
                Q_ASSERT(kk.size()==w.size());
            }
        }
#endif
        (*this)[w] = a + n;
        if(a+n == 0) {
            //Map::remove(w);
        }
        check();
        return a+n;
    }

    int total() const {
        return total_;
    }
private:
    int total_;
    int order_;
};



#endif // QCOUNTER_H
