#ifndef QQUTIL_H
#define QQUTIL_H


inline QTextStream& qqStdOut() {
    static QTextStream ts( stdout );
    return ts;
}


inline QStringList qqGlob(const QString& pattern) {
    //DEBUG("qGlob " << pattern)
    QFileInfo fi(pattern);
    QDirIterator it(fi.dir().absolutePath(), QStringList() << fi.fileName(), QDir::Files);
    QStringList filepaths;
    while (it.hasNext()) {
        filepaths.push_back(it.next());
    }
    return filepaths;
}


template<typename D>
QString dumps(const D& vec, QString sep=",") {
    QString result = "[";
    QTextStream ts(&result);
    for(auto it = vec.begin(); it!=vec.end(); it++) {
        if(it!=vec.begin())
            ts << sep;
        ts << *it;
    }
    ts<<"]";
    return result;
}

template<typename D,
         typename decltype(D::mapped_type)::type = 0>
QString dumps(const D& dict, QString sep=",") {
    QString result = "{";
    QTextStream ts(&result);
    for(auto it: QQItems<D>(dict)) {
        if(result.length()>1)
            ts << sep;
        ts << " " << it.key() << ": " << it.value();
    }
    ts<<"}";
    return result;
}

#endif // QQUTIL_H
