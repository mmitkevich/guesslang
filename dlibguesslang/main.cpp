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



int main()
{
    typedef char char_t;
    utf8_uifstream in("/dev/stdin");

    std::cout << "Hello\n";
    std::locale loc("en_US.UTF-8");

    uchar c = u'Я';
    bool a = utf8::isalpha(c);
    std::cout << format("a=%d %d %d",a, c, sizeof(c))<<"\n";
    while(1) {
        ustring s = read_utf8_n(in, 2, utf8::isalpha);
        auto u8str = convert_utf32_to_utf8(s);
        std::cout << format("\nREAD:[%s]\n", u8str.c_str());
    }
    return 0;
}
