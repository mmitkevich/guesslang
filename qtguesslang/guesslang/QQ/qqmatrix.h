#ifndef QQMATRIX_H
#define QQMATRIX_H

#include <algorithm>

template<typename T>
class QQVectorMatrixData: public QVector<T>, public QSharedData
{
public:
    typedef T value_type;
    using QVector<T>::QVector;
    using QVector<T>::begin;
    using QVector<T>::end;

    QQVectorMatrixData(int nrow, int ncol, T* d=nullptr) :
        QVector<T>(nrow*ncol)
    {
        nrow_ = nrow;
        ncol_ = ncol;
        QVector<T>::reserve(nrow_*ncol_);
        int i;
        if(d)
            for(i=0; i<ncol*nrow; i++)
                data()[i] = d[i];
        else
            for(i=0; i<ncol*nrow; i++)
                data()[i] = T();
    }

    int nr() const {
        return nrow_;
    }

    int nc() const {
        return ncol_;
    }

    T& operator()(int irow, int icol) {
        return (*this)[icol*nrow_+irow];
    }

    const T& operator()(int irow, int icol) const {
        return (*this)[icol*nrow_+irow];
    }

    T sumc(int irow, int ibegincol=0, int iendcol=-1) const {
        T result = T();
        ibegincol = ibegincol % nc();
        iendcol = iendcol % nc();
        if(ibegincol>iendcol)
            std::swap(ibegincol, iendcol);
        for(auto it = begin()+ibegincol*nrow_+irow; it <= begin()+iendcol*nrow_+irow; it += nrow_) {
            result += *it;
        }
        return result;
    }

    T sumr(int icol, int ibeginrow=0, int iendrow=-1) const {
        T result = T();
        ibeginrow = ibeginrow % nr();
        iendrow = iendrow % nr();
        if(ibeginrow>iendrow)
            std::swap(ibeginrow, iendrow);
        for(auto it = begin()+icol*nrow_+ibeginrow; it <= begin()+icol*nrow_+iendrow; it += 1) {
            result += *it;
        }
        return result;
    }

public:
    int nrow_;
    int ncol_;
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
        QQMatrix(new D(nrow, ncol), nrow, ncol, 0, 0) {
    }


    QQMatrix(int nrow, int ncol, const T* data) :
        QQMatrix(new D(nrow, ncol, data), nrow, ncol, 0, 0) {
    }

    QQMatrix(D* data, int nrow, int ncol, int irow, int icol):
        d(data),
        nrow_(nrow),
        ncol_(ncol),
        irow_(irow),
        icol_(icol) {
    }

    typedef typename D::iterator iterator;
    typedef typename D::const_iterator const_iterator;
    typedef typename D::value_type mapped_type;
    typedef int key_type;

    int nr() const {
        return nrow_;
    }

    int nc() const {
       return ncol_;
    }

    QQMatrix<T,D> subMatrix(int irowbegin, int icolbegin, int irowend, int icolend) {
        return QQMatrix(d, irowend-irowbegin+1, icolend-icolbegin+1, irow_+irowbegin, icol_+icolbegin);
    }

    QQMatrix<T,D> rows(int irowbegin, int irowend) {
       return QQMatrix(d, irowend-irowbegin+1, nr(), irow_+irowbegin, icol_);
    }

    QQMatrix<T,D> columns(int icolbegin, int icolend) {
      return QQMatrix(d, nr(), icolend-icolbegin+1, irow_, icol_+icolbegin);
    }

    QQMatrix<T,D> sumr() {
        QQMatrix<T,D> result(1, ncol_, 0, 0);
        for(int icol=0; icol<ncol_; icol++) {
            result(0, icol) += d->sumr(icol_+icol,irow_, irow_+nrow_-1);
        }
        return result;
    }

    QQMatrix<T,D> sumc() {
        QQMatrix<T,D> result(nrow_, 1, 0, 0);
        for(int irow=0; irow<nrow_; irow++) {
            result(irow, 0) += d->sumc(irow_+irow, icol_, icol_+ncol_-1);
        }
        return result;
    }


    T& operator()(int irow, int icol) {
        return (*d)(irow, icol);
    }

    const T& operator()(int irow, int icol) const {
        return (*d)(irow, icol);
    }

protected:
    QSharedDataPointer<D> d;
    int nrow_;
    int irow_;
    int ncol_;
    int icol_;
};


#endif // QQMATRIX_H
