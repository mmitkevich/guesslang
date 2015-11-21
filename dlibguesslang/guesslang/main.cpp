#include "deps.h"
#include "guesslang.h"
#include <fstream>
#include "utf8pp.h"

using namespace dlib;
using namespace guesslang;

std::function<bool(int)> regex_filter(ustring include)
{
    //std::regex re(include);
    //return [re](int ch){ return regex_match(ch, )
    return nullptr;
}



int test_read_n()
{
    typedef char char_t;
    utf8_uifstream in("/dev/stdin");

    std::cout << "Hello\n";
    std::locale loc("en_US.UTF-8");

    uchar c = u'Я';
    bool a = utf8::isalpha(c);
    std::cout << utf8::format("a=%d %d %d",a, c, sizeof(c))<<"\n";
    counter cnt;
    while(1) {
        ustring s = read_n(in, 2, utf8::isalpha);
        auto u8str = convert_utf32_to_utf8(s);
        std::cout << utf8::format("\nREAD:[%s]\n", u8str.c_str());
        cnt.push(s);
        for(auto kv: cnt){
            std::cout << convert_utf32_to_utf8(kv.first).c_str() << "=" << kv.second;
        }
    }
    return 0;
}

int main()
{
    guesslang::counter counter(2);
    dlib::utf8_uifstream ucin("/dev/stdin");
    counter.read_all(ucin);
    return 0;
}
