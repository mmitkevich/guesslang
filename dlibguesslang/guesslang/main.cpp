#include "deps.h"
#include "guesslang.h"
#include <fstream>


using namespace dlib;
using namespace guesslang;

std::function<bool(int)> regex_filter(ustring include)
{
    //std::regex re(include);
    //return [re](int ch){ return regex_match(ch, )
    return nullptr;
}





QTextStream & operator<<(QTextStream &stream, QStringList &strList) {
    stream << "[";
    bool first = true;
    for(auto s:strList){
        stream << s;
        if(first){
            stream << ", ";
            first = false;
        }
    }
    stream << "]";
    return stream;
}

int main(int argc, char* argv[])
{
    qStdOut() << "hello" << endl;

    struct options {
        QString training_glob;
        QString validation_glob;
        int k_max;
    };

    options opts;

    QCommandLineParser parser;
    parser.setApplicationDescription("Guess Language");
    parser.addOption(QCommandLineOption("l", "learn files mask", "learn", "*.sample"));
    parser.addOption(QCommandLineOption("v", "validate files mask", "validate", "*.sample"));



    QCoreApplication app(argc, argv);
    try {
        parser.process(app);
    }catch(std::exception &e) {
        qStdOut()<<"failed to parse command line:" << e.what();
        return -1;
    }

#ifndef NDEBUG
    for(auto opt: parser.optionNames()) {
        qStdOut() << opt << " = " << parser.value(opt) << endl;
    }
#endif

    QStringList fns, validate_fns;

    DEBUG("cwd="<<QDir::currentPath())

    if(parser.isSet("l")) {
        QString val = parser.value("l");
        fns << qGlob(val);
        DEBUG("found "<< fns.size() << " learn samples in "<<val);
    }
    if(parser.isSet("v")) {
        QString val = parser.value("v");
        validate_fns << qGlob(val);
        DEBUG("found "<< validate_fns.size() << " validate samples in "<<val);
        fns << validate_fns;
    }

    const int klassifier_order = 2;

    guesslang::QClassifier klassifier(2, 15);
    for(const auto &fn: fns) {
        QFile file(fn);
        file.open(QIODevice::ReadOnly);
        QTextStream is(&file);
        qStdOut() << "reading "<<fn<<endl;
        QFileInfo fi(fn);
        int n_chars;
        do {
            QParagraph sample(klassifier_order, guesslang::is_alpha, fi.baseName());
            n_chars = sample.read(is,"",100,1000);
            if(n_chars>0)
            {
                klassifier.append(sample);
                int n_itrs = klassifier.shuffle(50);
                DEBUG("S"<<klassifier.n_samples()<<" | " << fi.baseName() << " | "<<n_chars<< " chr | " << sample.str()<<" | LL="<<klassifier.likelihood()<<" | n_lang "<<klassifier.size()<<" | n_learn "<<fns.size()-validate_fns.size()<<" | n_validate "<<validate_fns.size());
                DEBUG(klassifier.str());

            }
        }while(n_chars>0);
    }
    //guesslang::QCounter counter(2, guesslang::is_alpha);

    //dlib::utf8_uifstream ucin("/dev/stdin");
    //QString s = QString::fromUtf8("abcdefg МИЩА");
    //QTextStream in(&s);
    //counter.read_all(in);
    return 0;
}
