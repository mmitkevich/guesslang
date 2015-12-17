//
// Created by mike on 11/20/15.
//

#ifndef GUESSLANG_CONFIG_H
#define GUESSLANG_CONFIG_H

#ifndef NDEBUG2
#define DEBUG(x)  { qStdOut() << x << endl; }
#else
#define DEBUG(x)
#endif

#define IF_DEBUG(x) x

#endif //GUESSLANG_CONFIG_H
