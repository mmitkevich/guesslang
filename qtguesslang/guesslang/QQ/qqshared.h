#ifndef QQSHARED_H
#define QQSHARED_H

#include <QtCore/QSharedData>
#include <QtCore/QSharedDataPointer>

template<typename D>
class QQShared {
public:
    QQShared(D* data):
        d(data) { }
    QQShared():
        d(new D){ }

    void swap(QQShared<D> &other) {
        std::swap(d, other.d);
    }

protected:
    D* data() { return d; }
    const D* data() const { return d; }
    QSharedDataPointer<D> d;
};

#endif // QQSHARED_H
