#ifndef GUESSLANG_CHARSTAT_H__
#define GUESSLANG_CHARSTAT_H__

#include "deps.h"
#include <regex>

namespace guesslang
{
    typedef dlib::ustring ustring;
    typedef utf8::uchar uchar;

    template<typename K, typename V>
    inline V value_of(std::map<K,V> dict, K key, V val=V()){
        auto it = dict.find(key);
        if (it==dict.end())
            return val;
        return it->second;
    }

    inline std::string convert_utf32_to_utf8(const ustring& u32str) {
        std::string u8str;
        for(auto it = u32str.begin(); it!=u32str.end(); ++it) {
            utf8::append (*it, std::back_insert_iterator<std::string>(u8str));
        }
        return u8str;
    }

    inline dlib::ustring read_n(dlib::utf8_uifstream &stream, int count,
                                std::function<bool(int)> filter = [](int ch) { return true; }, uchar pad = ' ') {
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
        typedef std::map<ustring, int> dict;

        /// counts character sequence of length m in input
        counter(int m = 1, std::function<bool(int)> filter = [](int ch) { return true; }):
            total_(0),
            filter_(filter),
            m_(m)
        {
            assert(m>0);
        }

        /// probability of charcter
        double probability(ustring c) {
//            auto it = cnt_.find(c);
//            return it==cnt_.end() ? 0. : double(it->second) / total_;
            return double(value_of(cnt_, c, 0))/total_;
        }

        inline void read_all(dlib::utf8_uifstream& stream) {
            uchar c;
            while((c = stream.get())!=EOF) {
                if(filter_(c))
                    put(c);
#ifndef NDEBUG
                else
                    std::cout<<"IGNORED "<<c<<"\n";
#endif
            }
        }

        counter& put(uchar c) {
            if(c!=EOF) {
                buf_.push_back(c);
                if (buf_.size() >= m_) {
                    ustring str(buf_.begin(), buf_.end());
                    push(str);
                }
            }
            return *this;
        }

        /// update counter
        void push(const dlib::ustring& s, int n=1) {
#ifndef NDEBUG
            std::cout<<"PUSH "<<s.c_str()<<"\n";
#endif
            total_ += n;
            cnt_[s] += n;
        }

        int total() const {
            return total_;
        }

        dict::const_iterator cbegin() const {
            return cnt_.cbegin();
        }

        dict::const_iterator cend() const {
            return cnt_.cend();
        }

        dict::iterator begin() {
            return cnt_.begin();
        }

        dict::iterator end() {
            return cnt_.end();
        }

        typedef dict::iterator iterator;

    private:
        std::deque<uchar> buf_;
        int m_;
        dict cnt_;
        int total_;
        std::function<bool(int)> filter_;
    };
};

#endif //GUESSLANG_CHARSTAT_H__
