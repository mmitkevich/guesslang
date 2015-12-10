#ifndef QQMATRIX_H
#define QQMATRIX_H

#include <algorithm>
#include <functional>
#include <limits>
#include <cstddef>

struct QQIndex2 {
    int c;
    int r;

    static const QQIndex2& null() {
        static QQIndex2 n(-1000000000,-1000000000);
        return n;
    }

    QQIndex2(int r, int c) {
        this->c = c;
        this->r = r;
    }

    int length() const{
        return c*r;
    }

    int trace() const {
        return c+r;
    }

    QQIndex2& operator+=(const QQIndex2 &other) {
        c+=other.c;
        r+=other.r;
        return *this;
    }

    QQIndex2 operator+(const QQIndex2 &other) const {
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

    static QQIndex2 axis(int axis, int len=1) {
        switch(axis){
        case 0:
            return QQIndex2(len, 0);
        default:
            return QQIndex2(0, len);
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

    QQIndex2 shape() const {
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
    typedef QVector<T> Vector;
    using Vector::QVector;
    using Vector::begin;
    using Vector::end;
    using Vector::data;
    using Vector::resize;
    using Vector::size;

    QQVectorMatrixData(const QQIndex2 &sz, const T* d=nullptr) :
        QVector<T>(sz.length()),
        sz_(sz)
    {
        resize(sz_.r, sz_.c);
        int i;
        int n = sz_.length();
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

    int offset(const QQIndex2 &ix) const {
        Q_ASSERT(sz_.contains(ix));
        int ofs = ix.offsetr(sz_.c);
        Q_ASSERT(ofs>=0);
        int s = size();
        Q_ASSERT(ofs<s);
        return ofs;
    }

    void ensure(int ir=-1, int ic=-1) {
        if(ic>=0)
            Q_ASSERT(ic<sz_.c);
        if(ir>=sz_.r){
            resize(ir+1, 0);
            QVector<T>::resize(sz_.length());
        }
    }

    void resize(int nr=0, int nc=0) {
        if(nc>0)
            Q_ASSERT(nc==sz_.c);

        if(nr>0) {
            sz_.r = nr;
            QVector<T>::resize(sz_.length());
        }

    }

    T& at(const QQIndex2 &ix) {
        ensure(ix.r, ix.c);
        T &v = (*this)[offset(ix)];
        return v;
    }

    const T& at(const QQIndex2 &ix) const {
        //ensure(ix.r, ix.c);
        const T& v = Vector::at(offset(ix));
        return v;
    }

    void set(const QQIndex2 &ix, const T&value) {
        int ofs = offset(ix);
        (*this)[ofs] = value;
        auto v2 = Vector::at(ofs);
    }

    template<typename F>
    void apply(const QQVectorMatrixData<T> &other, const QQRange2& range, F& func) {
        for(QQIndex2 ix=range.lt; ix.r<=range.rb.r; ix.r++)
            for(ix.c=range.lt.c; ix.c<=range.rb.c; ix.c++)
                at(ix) = func(at(ix), other.at(ix));
    }

    template<typename F, typename R>
    R reduce(const QQRange2& range, F& func, R start=R()) const {
        for(QQIndex2 ix=range.lt; ix.r<=range.rb.r; ix.r++)
            for(ix.c=range.lt.c; ix.c<=range.rb.c; ix.c++)
                start = func(at(ix), start);
        return start;
    }

    template<typename F, typename R>
    QQIndex2 argmax(const QQRange2& range, F& func, R start=R()) const {
        QQIndex2 ix;
        QQIndex2 jx(-1, -1);
        for(QQIndex2 ix=range.lt; ix.r<=range.rb.r; ix.r++)
            for(ix.c=range.lt.c; ix.c<=range.rb.c; ix.c++) {
                R val = func(at(ix));
                if(val>=start){
                    jx = ix;
                    start = val;
                }
            }
        return jx;
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
template<typename T=double, typename D = QQVectorMatrixData<T> >
class QQMatrix
{
public:
    QQMatrix(int nrow, int ncol):
        QQMatrix(QQIndex2(nrow, ncol), nullptr) {
    }


    QQMatrix(const QQIndex2 &sz, const T* data) :
        QQMatrix(new D(sz, data), sz) {
    }

    QQMatrix(int nrow, int ncol, const std::initializer_list<T>& list) :
        QQMatrix(new D(QQIndex2(nrow,ncol), list.begin()), QQIndex2(nrow,ncol)) {
    }

    QQMatrix(D* data, const QQRange2 &range):
        d(data),
        range_(range) {
    }

    QQMatrix(QQMatrix&& other): range_(std::move(other.range_))
    {
        d.swap(other.d);
    }

    QQMatrix(const QQMatrix& other):
        d(other.d),
        range_(other.range_){

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

    QQIndex2 shape() const {
        return range_.shape();
    }

    QQMatrix& operator=(const QQMatrix &other) {
        Q_ASSERT(nr()==other.nr());
        Q_ASSERT(nc()==other.nc());
        for(int ir=0;ir<other.nr();ir++)
            for(int ic=0;ic<other.nc();ic++)
                at(ir,ic) = other(ir,ic);
        return *this;
    }

    QQMatrix& operator=(QQMatrix &&other) {
        //d.swap(other.d);
        //std::swap(range_, other.range_);
        Q_ASSERT(nr()==other.nr());
        Q_ASSERT(nc()==other.nc());
        for(int ir=0;ir<other.nr();ir++)
            for(int ic=0;ic<other.nc();ic++)
                at(ir,ic) = other(ir,ic);
        return *this;
    }

    bool operator==(const QQMatrix<T,D> &other) const {
        Q_ASSERT(nr()==other.nr());
        Q_ASSERT(nc()==other.nc());
        for(int ir=0;ir<other.nr();ir++)
            for(int ic=0;ic<other.nc();ic++)
            {
                auto v1 = at(ir, ic);
                auto v2 = other.at(ir,ic);
                if(v1 != v2)
                    return false;
            }
        return true;
    }

    QQMatrix<T,D> subm(const QQRange2 &range) {
        return QQMatrix(d, range_.subrange(range));
    }

    QQMatrix<T,D> &setsubm(const QQIndex2 &lt, const QQMatrix<T,D>& other) {
        for(int ir=0; ir<other.nr(); ir++)
            for(int ic=0; ic<other.nc(); ic++)
                at(ir+lt.r, ic+lt.c) = other(ir,ic);
        return *this;
    }

    QQRange2 axis(int axis, int i) {
        auto range = range_.axis(axis);
        range += QQIndex2::axis(!axis, i);
        return range;
    }

    void resize(int nr, int nc=-1) {
        d->resize(nr, nc);
    }

    QQMatrix<T,D> row(int i) {
        QQRange2 range = axis(1, i);
        return QQMatrix<T,D>(d, range);
    }

    QQMatrix<T,D> col(int i) {
      return QQMatrix<T,D>(d, axis(0, i));
    }

    QQMatrix<T,D>& setrow(int i, const QQMatrix<T,D> &other) {
        setsubm(QQIndex2(i,0), other);
        return *this;
    }

    QQMatrix<T,D>& setcol(int i, const QQMatrix<T,D> &other) {
        setsubm(QQIndex2(0,i), other);
        return *this;
    }

    template<typename F, typename R>
    R reduce(const R& start, const F& func) const {
        return d->reduce(range_, func, start);
    }

    T max() const {
        return reduce(std::numeric_limits<T>::min(),
            [&](const T& item, const T& start)->T {
                return std::max(item, start);
            });
    }

    template<typename F, typename R>
    QQIndex2 argmax(R& start=R(), F& func = F()) const {
        return d->argmax(range_, func, start);
    }

    const T& at(int irow, int icol) const {
        return at(QQIndex2(irow, icol));
    }

    T& at(int irow, int icol) {
        return at(QQIndex2(irow, icol));
    }

    T& at(const QQIndex2& ix) {
        T& value = d->at(range_.lt+ix);
        // fixup the range
        ensure(ix);
        return value;
    }

    const T& at(const QQIndex2& ix) const {
        const T& value = d->at(range_.lt+ix);
        return value;
    }

    void set(int irow, int icol, const T&value)  {
        d->set(range_.lt+QQIndex2(irow, icol), value);
    }

    void ensure(const QQIndex2& ix){
        if(ix.r>=range_.rb.r)
            range_.rb.r = ix.r;
        if(ix.c>=range_.rb.c)
            range_.rb.c = ix.c;
    }

    template<typename F>
    QQMatrix<T,D> reduce_axis(int axis, F&func) const {
        QQRange2 subrange = range_.axis(1-axis);
        QQMatrix<T,D> result(subm(subrange));
        while(range_.contains(subrange+=QQIndex2::axis(axis))){
            func(result, subm(subrange));
        }
        return result;
    }

    QQMatrix<T,D> sum_axis(int axis) const {
        return binary(axis,[&](QQMatrix<T,D> &result, const QQMatrix<T,D> &other) { result += other; });
    }

    QQMatrix<T,D>& operator+=(const QQMatrix<T,D> &other) {
        Q_ASSERT(shape() == other.shape());
        d->apply(other->d, range_, [&](const T&a, const T&b){return a+b;});
        return *this;
    }

    QQMatrix<T,D> operator+(const QQMatrix<T,D> &other) {
        Q_ASSERT(shape() == other.shape());
        QQMatrix<T,D> result(*this);
        result += other;
        return result;
    }

    T& operator()(int irow, int icol) {
        return at(QQIndex2(irow, icol));
    }

    const T& operator()(int irow, int icol) const {
        return at(QQIndex2(irow, icol));
    }

    T& operator[](int ioffset) {
        return (*d)[d->offset(range_.lt) + ioffset];
    }

protected:
    QSharedDataPointer<D> d;
    QQRange2 range_;
};

template<int NC, typename T = double, typename D = QQVectorMatrixData<T>>
static QQMatrix<T,D> qqmat(const std::initializer_list<T>& arr) {
    return QQMatrix<T,D>(QQIndex2(arr.size()/NC + ((arr.size()%NC)!=0), NC), arr.begin());
}

template<typename T, typename D>
QTextStream& operator<<(QTextStream &s, const QQMatrix<T,D> &mat) {
    for(int ir=0; ir<mat.nr(); ir++) {
        for(int ic=0; ic<mat.nc(); ic++) {
            if(ic>0)
                s << "\t";
            s << mat.at(ir, ic);
        }
        s << "\n";
    }
    return s;
}

#endif // QQMATRIX_H
