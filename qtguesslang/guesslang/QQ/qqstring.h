#ifndef QQSTRING_H
#define QQSTRING_H

class QTextStream;
class QStringList;

QTextStream & operator<<(QTextStream &stream, QStringList &strList);

/// QString outstr;
/// outstr << QList<int>({1,2});
/// outstr == "[1,2]"
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
    stream << "]";
    return stream;
}

/// QString outstr = "one=";
/// outstr << 1;
/// outstr == "one=1"
template<typename T>
inline QString& operator<<(QString &str, const T& arg) {
    str.append(arg);
    return str;
}


/// QString("%1=%2") % "123" % 123;
template<typename T>
inline QString operator%(const QString &str, const T &arg) {
    return str.arg(arg);
}

#endif // QQSTRING_H
