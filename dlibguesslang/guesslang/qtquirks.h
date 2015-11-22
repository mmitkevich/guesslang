#ifndef QTQUIRKS_H
#define QTQUIRKS_H

QTextStream& qStdOut();

QStringList qGlob(const QString& pattern);


template<typename T>
QTextStream& operator<<(QTextStream& stream, const QList<T> &list) {
    stream << "[";
    bool first = true;
    for(int i=0; i<list.size(); i++) {
        if(!first)
            stream << ", ";
        stream << list[i];
        first = false;
    }
    return stream;
}

#endif // QTQUIRKS_H
