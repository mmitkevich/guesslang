#ifndef QQMATRIX_H
#define QQMATRIX_H

#include <algorithm>

struct QQIndex2 {
    int c;
    int r;

    QQIndex2(int r, int c) {
        this->c = c;
        this->r = r;
    }

    int length() const{
        return c*r;
    }

    QQIndex2& operator+=(const QQIndex2 &other) {
        c+=other.c;
        r+=other.r;
        return *this;
    }

    QQIndex2 operator+(const QQIndex2 &other) {
        QQIndex2 result(*this);
        result+=other;
        return result;
    }

    int offsetr(int ncol) const {
        return c+ncol*r;
    }

    int offsetc(int nrow) const {
        return r+nrow*c;
    }

    bool contains(const QQIndex2 &ix) const {
        return ix.c>=0 && ix.c<c && ix.r>=0 && ix.r<r ;
    }

    QQIndex2 lt(const QQIndex2 &other) const {
        return QQIndex2(std::min(r,other.r),std::min(c,other.c));
    }

    QQIndex2 rb(const QQIndex2 &other) const {
        return QQIndex2(std::max(r,other.r),std::max(c,other.c));
    }

    QQIndex2 operator%(const QQIndex2& sz) const {
        return QQIndex2(r%sz.r, c%sz.c);
    }

    static QQIndex2 axis(int axis) {
        switch(axis){
        case 0:
            return QQIndex2(1, 0);
        default:
            return QQIndex2(0, 1);
        }
    }
};

struct QQRange2 {
    QQIndex2 lt;
    QQIndex2 rb;

    QQRange2(int r1, int c1, int r2, int c2):
        lt(r1,c1),
        rb(r2,c2){
    }

    QQRange2(const QQIndex2 &ix1, const QQIndex2 &ix2):
        lt(ix1),
        rb(ix2){
    }

    QQRange2(const QQIndex2 &sz):
        lt(0, 0),
        rb(sz.r-1, sz.c-1) {
    }

    int length(int axis) const {
        switch(axis){
        case 0: return rb.r-lt.r+1;
        case 1: return rb.c-lt.c+1;
        }
        throw std::bad_exception();
    }

    QQIndex2 size() const {
        return QQIndex2(rb.r-lt.r+1, rb.c-lt.c+1);
    }

    QQRange2 axis(int axis) const{
        switch(axis){
        case 0:
            return QQRange2(lt.r, lt.c, rb.r, lt.c);
        default:
            return QQRange2(lt.r, lt.c, lt.r, rb.c);
        }
    }

    QQRange2 subrange(const QQRange2& inner) {
        return QQRange2(inner.lt.rb(lt).lt(rb), inner.rb.lt(rb).rb(lt));
    }

    bool contains(const QQIndex2 &ix) {
        return ix.c>=lt.c && ix.c<=rb.c && ix.r>=lt.r && ix.r<=rb.r;
    }

    QQRange2& operator+=(const QQIndex2 &other) {
        lt+=other;
        rb+=other;
        return *this;
    }

    QQRange2 operator+(const QQIndex2 &other) {
        QQRange2 result(*this);
        result += other;
        return result;
    }

};

template<typename T>
class QQVectorMatrixData: public QVector<T>, public QSharedData
{
public:
    typedef T value_type;
    using QVector<T>::QVector;
    using QVector<T>::begin;
    using QVector<T>::end;

    QQVectorMatrixData(const QQIndex2 &sz, T* d=nullptr) :
        QVector<T>(sz.length()),
        sz_(sz)
    {
        int n = sz_.length();
        QVector<T>::reserve(n);
        int i;
        if(d)
            for(i=0; i<n; i++)
                QVector<T>::data()[i] = d[i];
        else
            for(i=0; i<n; i++)
                QVector<T>::data()[i] = T();
    }

    int nr() const {
        return sz_.r;
    }

    int nc() const {
        return sz_.c;
    }

    int offset(const QQIndex2 &ix) {
        Q_ASSERT(sz_.contains(ix));
        return ix.offsetr(sz_.c);
    }

    void reserve(int nrow) {
        if(nrow>sz_.r){
            sz_.r = nrow+1;
            QVector<T>::reserve(sz_.length());
        }
    }

    T& at(const QQIndex2 &ix) {
        reserve(ix.r+1);
        return (*this)[offset(ix)];
    }

    template<typename F>
    void apply(const QQVectorMatrixData<T> &other, F&func) {
        QQIndex2 ix;
        QQRange2 range(sz_);
        for(QQIndex2 ix=range.lt; ix.r<=range.rb.r; ix.r++)
            for(ix.c=range.lt.c; ix.c<=range.rb.c; ix.c++)
                at(ix) = func(at(ix), other.at(ix));
    }

    template<typename F, typename R>
    R reduce(const QQRange2& range, F& func, R start=R()) const {

        QQIndex2 ix;
        for(QQIndex2 ix=range.lt; ix.r<=range.rb.r; ix.r++)
            for(ix.c=range.lt.c; ix.c<=range.rb.c; ix.c++)
                start = func(at(ix), start);
        return start;
    }


public:
    QQIndex2 sz_;
};
#if 0

template<typename T>
class QQListVectorMatrixData: public QList<QVector<T> >, public QSharedData
{
public:
    typedef T value_type;
    using QList<QVector<T> >::QList;
    using QList<QVector<T> >::begin;
    using QList<QVector<T> >::end;

    QQListVectorMatrixData(int nrow, int ncol, const T* data=nullptr)
    {
        QList<QVector<T>>::reserve(ncol);
        int irow, icol;
        for(icol=0; icol<ncol; icol++) {
            QVector<T> &vec = (*this)[icol];
            for(irow=0; irow<nrow; irow++)
                vec[irow] = data ? data[ncol*irow+icol] : T();
    }

    int nr() const {
        return (*this)[0].size();
    }

    int nc() const {
        return (*this).size();
    }

    T& operator()(int irow, int icol) {
        return (*this)[icol][irow];
    }

    const T& operator()(int irow, int icol) const {
        return (*this)[icol][irow];
    }

    T sumc(int irow, int ibegincol=0, int iendcol=-1) const {
        T result = T();
        ibegincol = ibegincol % nc();
        iendcol = iendcol % nc();
        if(ibegincol>iendcol)
            std::swap(ibegincol, iendcol);
        for(int icol = ibegincol; icol <= iendcol; icol++) {
            result += (*this)[icol][irow];
        }
        return result;
    }

    T sumr(int icol, int ibeginrow=0, int iendrow=-1) const {
        T result = T();
        ibeginrow = ibeginrow % nr();
        iendrow = iendrow % nr();
        if(ibeginrow>iendrow)
            std::swap(ibeginrow, iendrow);
        QVector<T> &vec = (*this)[icol];
        for(int irow = ibeginrow; irow <= iendrow; irow++) {
            result += vec[irow];
        }
        return result;
    }
};
#endif
/// submatrix
template<typename T, typename D>
class QQMatrix
{
public:
    QQMatrix(int nrow, int ncol):
        QQMatrix(QQIndex2(nrow, ncol), nullptr) {
    }


    QQMatrix(const QQIndex2 &sz, const T* data) :
        QQMatrix(new D(sz, data), sz) {
    }

    QQMatrix(D* data, const QQRange2 &range):
        d(data),
        range_(range) {
    }

    typedef typename D::iterator iterator;
    typedef typename D::const_iterator const_iterator;
    typedef typename D::value_type mapped_type;
    typedef int key_type;

    int nr() const {
        return range_.length(0);
    }

    int nc() const {
       return range_.length(1);
    }

    QQMatrix<T,D> subm(const QQRange2 &range) {
        return QQMatrix(d, range_.subrange(range));
    }

    QQMatrix<T,D> axis(int axis, int i) {
        return QQMatrix(d, range_.axis(axis)+QQIndex2::axis(!axis));
    }

    QQMatrix<T,D> row(int i) {
       return axis(0, i);
    }

    QQMatrix<T,D> col(int i) {
      return axis(1, i);
    }

    template<typename F, typename R>
    R reduce(F& func, R start=R()) const {
        return d->reduce(range_, func, start);
    }

    T max(const QQRange2 &r) const {
        return reduce(range_, [&](const T& item, const T& start){return std::max(item, start);});
    }

    T& at(int irow, int icol) {
        return at(QQIndex2(irow, icol));
    }

    T& at(const QQIndex2& ix) {
        return d->at(ix);
    }

    template<typename F>
    QQMatrix<T,D> reduce(int axis, F&func) const {
        QQRange2 subrange = range_.axis(1-axis);
        QQMatrix<T,D> result(subm(subrange));
        while(range_.contains(subrange+=QQIndex2::axis(axis))){
            func(result, subm(subrange));
        }
        return result;
    }

    QQMatrix<T,D> sum(int axis) const {
        return binary(axis,[&](QQMatrix<T,D> &result, const QQMatrix<T,D> &other) { result += other; });
    }

    QQMatrix<T,D>& operator+=(const QQMatrix<T,D> &other) {
        d->apply(other->d, [&](const T&a, const T&b){return a+b;});
        return *this;
    }

    QQMatrix<T,D> operator+(const QQMatrix<T,D> &other) {
        QQMatrix<T,D> result(*this);
        result += other;
        return result;
    }

    T& operator()(int irow, int icol) {
        return d->at(irow, icol);
    }

    const T& operator()(int irow, int icol) const {
        return d->at(irow, icol);
    }

protected:
    QSharedDataPointer<D> d;
    QQRange2 range_;
};


#endif // QQMATRIX_H
