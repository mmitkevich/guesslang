#include "deps.h"
#include "guesslang.h"
#include <fstream>


int main(int argc, char* argv[])
{
    qqStdOut() << "hello" << endl;

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
    parser.addOption(std::move(QCommandLineOption("l", "*.learn", "*.learn")));
    parser.addOption(std::move(QCommandLineOption("v", "*.validate", "*.validate")));



    QCoreApplication app(argc, argv);
    parser.process(app);

#ifndef NDEBUG
    for(auto opt: parser.optionNames()) {
        qqStdOut() << opt << " = " << parser.value(opt) << endl;
    }
#endif

    QStringList learn_fns, validate_fns;

    DEBUG("cwd="<<QDir::currentPath())

    if(parser.isSet("l")) {
        learn_fns << qqGlob(parser.value("l"));
        DEBUG("found "<< learn_fns.size() << " learn samples in "<<parser.value("l"));
    }
    if(parser.isSet("v")) {
        validate_fns << qqGlob(parser.value("v"));
        DEBUG("found "<< validate_fns.size() << " validate samples in "<<parser.value("v"));
    }

    learn_fns << validate_fns;

    const int klassifier_order = 2;
    const int max_lang = 15;
    const int max_shuffle = 10;
    const int nlsep = 1; // number of newlines in row to separate sampmles. 1=\n, 0=off, 2=\n\n

    QKMedoidsClassifier<QParagraph> klassifier(max_lang);
    for(auto fn: learn_fns) {
        QFile file(fn);
        QFileInfo fi(fn);
        file.open(QIODevice::ReadOnly);
        QTextStream is(&file);
        //qStdOut() << "reading " << fn << "\n";
        QVector<QParagraph> samples;
        int n_samples = QParagraph::read_all(is, samples, klassifier_order, QParagraph::is_alpha, 0, 1000, 10000);
        klassifier.append(samples);
        klassifier.shuffle();
        int n_itrs = 0;// klassifier.shuffle(max_shuffle);
        DEBUG(QString("S%1 | %2(%3) | ll=%4 | k=%5") % klassifier.samples().size() % fi.baseName() % n_samples
              % klassifier.likelihood() % klassifier.centroids().size());
        DEBUG(klassifier.str());
    }
    return 0;
}
