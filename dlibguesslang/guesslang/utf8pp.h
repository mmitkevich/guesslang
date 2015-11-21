//
// Created by mike on 11/19/15.
//

#ifndef GUESSLANG_UTF8PROC_H
#define GUESSLANG_UTF8PROC_H

#include <string>
#include <utf8proc.h>
#include <utf8.h>

namespace utf8 {

    inline bool isalpha(int c)
    {
        int cat = ::utf8proc_get_property(c)->category;
        return (UTF8PROC_CATEGORY_LU <= cat && cat <= UTF8PROC_CATEGORY_LO);
    }

// is_print    (c == 0x0601 || c == 0x0602 || c == 0x0603 || c == 0x06dd);

    typedef std::basic_string<int> ustring;
    typedef int uchar;

    template<typename ... Args>
    inline std::string format( const std::string& format, Args ... args )
    {
        size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
        std::unique_ptr<char[]> buf( new char[ size ] );
        snprintf( buf.get(), size, format.c_str(), args ... );
        return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
    }

}
#endif //GUESSLANG_UTF8PROC_H
