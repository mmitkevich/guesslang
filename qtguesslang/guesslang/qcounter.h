#ifndef QCOUNTER_H
#define QCOUNTER_H

#include <QtCore/QMap>

/// character probabilities
template<typename T>
class QCounter : public QMap<T, int>
{
public:
    //typedef QHash<QString, int> Dict;
    //typedef Dict::iterator iterator;
    //typedef Dict::const_iterator const_iterator;

    /// counts character sequences
    QCounter():
        total_(0)
    {
    }

    /// probability of charcter sequence
    double probability(const T &s) {
        return double(value(s, 0)) / total_;
    }

    /// add other counter into this counter
    QCounter& operator+=(const QCounter &other) {
        for(auto it: items(other)) {
            increment(it.key(), it.value());
        }
        return *this;
    }

    QCounter operator+(const QCounter &other) const {
        QCounter result = *this;
        result+=other;
        return result;
    }

    /// subtract rhs counter from this counter
    QCounter& operator-=(const QCounter &other) {
        for(auto it: items(other)) {
            increment(it.key(), -it.value());
        }
        return *this;
    }

    QCounter operator-(const QCounter &other) const {
        QCounter result = *this;
        result-=other;
        return result;
    }

    /// increment counter for sequence specified
    int increment(T& w, int n=1) {
        //DEBUG("inc " << s)
        total_ += n;
        (*this)[w] += n;
        auto v = (*this)[w];
        if(v==0) {
            QMap<QString,int>::remove(w);
        }
        return v;
    }

    int total() const {
        return total_;
    }
private:
    int total_;
};



#endif // QCOUNTER_H
