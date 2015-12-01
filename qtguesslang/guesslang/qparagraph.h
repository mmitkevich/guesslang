#ifndef QPARAGRAPH_H
#define QPARAGRAPH_H

#include <QtCore/QSharedData>
#include <QtCore/QSharedDataPointer>
#include <QtCore/QMap>

#include "qcounter.h"

typedef bool (*QCharPredicate)(const QChar&);
//template<typename V, typename T> V StaticCast(const T& d) {
//    return static_cast<V>(d);
//}

template<typename T>
bool True(const T& arg) { return true; }

class QParagraph : public QMap<int, QCounter<QString> >
{
public:
    typedef QMap<int, QCounter<QString>> base_type;

    QParagraph()
    {
    }

    static bool is_alpha(const QChar &c) {
        return QChar(c).isLetter();
    }

    int alphabet_size() const {
        return (*this)[1].size();
    }

    int order() const {
        return lastKey();
    }

    const QCounter<QString> at(int i) const {
        return (*this)[i];
    }

    QQItems<base_type> items() {
        return QQItems<base_type>(*static_cast<base_type*>(this));
    }

    QString str() const {
        QString r;
        return r << QString("par(")
                << (QString("%1, %2, %3") % at(2).total() % at(1).total() % at(1).size())
                << ")";
    }

    static int read_all(QTextStream &stream, QVector<QParagraph> &samples, int order, QCharPredicate filter = True<QChar>, int n_sep=1, int min_chars=0, int max_chars=10000, int *skipped=nullptr)
    {
        QQueue<QChar> buf;
        QString w;
        QParagraph sample;
        int n_samples = 0;
        int n_sepcnt = 0;
        int n_chars = 0;
        while(!stream.atEnd()) {
            QChar c;
            stream >> c;
            bool flush = n_chars>=max_chars;
            if(n_sep>0 && n_chars>=min_chars) {
                if(c=='\n')
                    n_sepcnt++;
                else
                    n_sepcnt = 0;
                if(n_sepcnt>=n_sep)
                    flush = true;
            }
            if(flush) {
                samples.append(std::move(sample));
                sample.clear();
                n_samples++;
                n_chars = 0;
            }
            if(filter(c))
            {
                buf.append(c.unicode());
                for(int i=1; i<=buf.size(); i++) {
                    QString w(i, ' ');
                    qCopy(buf.begin(), buf.end(), w.begin());
                    sample[i].increment(w);
                }

                if (buf.size() >= order) {
                    buf.pop_front();
                }
                n_chars++;
            }
            else if(skipped)
                *skipped++;
        }
        return n_samples;
    }
};

template<typename T>
inline QTextStream& operator<<(QTextStream& ts, const T& p) {
    ts << p.str();
    return ts;
}


#endif // QPARAGRAPH_H
