#ifndef GUESSLANG_CHARSTAT_H__
#define GUESSLANG_CHARSTAT_H__

#include "deps.h"
#include <regex>

namespace guesslang
{
    typedef dlib::ustring ustring;
    typedef dlib::unichar uchar;

    template<typename ... Args>
    inline std::string format( const std::string& format, Args ... args )
    {
        size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
        std::unique_ptr<char[]> buf( new char[ size ] );
        snprintf( buf.get(), size, format.c_str(), args ... );
        return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
    }

    inline std::string convert_utf32_to_utf8(const ustring& u32str) {
        std::string u8str;
        for(auto it = u32str.begin(); it!=u32str.end(); ++it) {
            utf8::append (*it, std::back_insert_iterator<std::string>(u8str));
        }
        return u8str;
    }

    inline dlib::ustring read_utf8_n(dlib::utf8_uifstream &stream, int count, std::function<bool(int)> filter = [](int ch){return true;}, uchar pad=' ') {
        dlib::ustring buf;
        buf.reserve(count);
        int l, c;
        while((l = buf.length()) < count
          && (c = stream.get())!=EOF) {
            if(filter(c))
                buf.push_back(c);
            else
                std::cout<<"IGNORED "<<c;
        }
        while(buf.length() < count) {
            buf.push_back(pad);
        }
        return buf;
    }

    /// character probabilities
    class counter
    {
    public:
        typedef std::map<dlib::ustring, int> dict;

        /// counts character sequence of length m in input
        counter(int m) {
            m_ = m;
        }

        /// update counter
        void read(dlib::utf8_uifstream &stream) {
            push(read_utf8_n(stream, m_));
        }
        void push(const dlib::ustring &s, int n=1) {
            cnt_[s]++;
        }
    private:
        dlib::ustring buf_;
        int m_;
        dict cnt_;
    };
};

#endif //GUESSLANG_CHARSTAT_H__
