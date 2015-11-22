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
    parser.setApplicationDescription("GuessLang");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addOption(std::move(QCommandLineOption(QStringList() << "learn" << "l", "*.learn", "*.learn")));
    parser.addOption(std::move(QCommandLineOption(QStringList() << "validate" << "v", "*.validate", "*.validate")));



    QCoreApplication app(argc, argv);
    parser.process(app);

#ifndef NDEBUG
    for(auto opt: parser.optionNames()) {
        qStdOut() << opt << " = " << parser.value(opt) << endl;
    }
#endif

    QStringList learn_fns, validate_fns;

    DEBUG("cwd="<<QDir::currentPath())

    if(parser.isSet("learn")) {
        learn_fns << qGlob(parser.value("learn"));
        DEBUG("found "<< learn_fns.size() << " learn samples in "<<parser.value("learn"));
    }
    //if(parser.isSet("validate"))
    //    validate_fns << qGlob(parser.value("validate"));

    const int klassifier_order = 2;

    guesslang::QClassifier klassifier(15);
    for(auto fn: learn_fns) {
        QFile file(fn);
        file.open(QIODevice::ReadOnly);
        QTextStream is(&file);
        qStdOut() << "reading "<<fn<<endl;
        QFileInfo fi(fn);
        int n_chars;
        do {
            QParagraph sample(klassifier_order, guesslang::is_alpha, fi.baseName());
            n_chars = sample.read(is,"",1000,10000);
            if(n_chars>0)
            {
                klassifier.append(sample);
                int n_itrs = klassifier.shuffle(10);
                DEBUG("S"<<klassifier.n_samples()<<",from " << fn << " read "<<n_chars<< " chars :" << sample.str()<<" fitness "<<klassifier.likelihood()<<" n_lang "<<klassifier.size());
                DEBUG("S"<<klassifier.n_samples()<<klassifier.str());

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
