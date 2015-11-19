#include <iostream>
#include "dlib/hash_map.h"
#include "dlib/any.h"

using namespace dlib;

typedef hash_map<std::string, int, 14>::kernel_1a Dict;

//template<typename key, typename value> using Dict = ;

void test_dlib()
{
    Dict  dict;

    std::string key = "A";
    int val = 1;
    dict.add(key, val);
    enumerable<map_pair<std::string,int>> &e = dict;
    enumerable<map_pair<std::string,int>> &e2 = dict;
    auto b1 = e.move_next();
    e2.reset();
    auto b2 =  e2.move_next();
    e2.reset();
    auto c1 =  e.move_next();
    //hash_map::kernel_1a;
    std::cout<<"HI";
}

//int main() {
//    test_dlib();
//    return 0;
//}
