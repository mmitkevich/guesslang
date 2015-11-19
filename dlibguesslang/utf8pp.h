//
// Created by mike on 11/19/15.
//

#ifndef GUESSLANG_UTF8PROC_H
#define GUESSLANG_UTF8PROC_H

#include <utf8proc.h>

namespace utf8 {

    inline bool isalpha(int c)
    {
        int cat = ::utf8proc_get_property(c)->category;
        return (UTF8PROC_CATEGORY_LU <= cat && cat <= UTF8PROC_CATEGORY_LO);
    }

// is_print    (c == 0x0601 || c == 0x0602 || c == 0x0603 || c == 0x06dd);

}
#endif //GUESSLANG_UTF8PROC_H
