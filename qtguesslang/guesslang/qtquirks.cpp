#include "deps.h"
#include "qtquirks.h"

QTextStream& qStdOut() {
    static QTextStream ts( stdout );
    return ts;
}


QStringList qGlob(const QString& pattern) {
    DEBUG("qGlob " << pattern);
    QFileInfo fi(pattern);
    QDirIterator it(fi.dir().absolutePath(), QStringList() << fi.fileName(), QDir::Files);
    QStringList filepaths;
    while (it.hasNext()) {
        filepaths.push_back(it.next());
    }
    return filepaths;
}

