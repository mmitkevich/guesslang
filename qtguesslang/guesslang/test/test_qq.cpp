#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../QQ/QQ.h"

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Matrix row should return correct row", "[QQMatrix]" ) {
    QQMatrix<double> m(2,2, {
        1., 2.,
        3., 4.});
    qqStdOut() << m;
    QQMatrix<double> m1(1,2, {1, 2});
    QQMatrix<double> m2 = m.row(0);
    qqStdOut() << m2;
    REQUIRE( m2 == m1 );
}

QQMatrix<double> return_matrix() {
    return qqmat<2>({1.,2.,3.,4.});
}

TEST_CASE("Can return matrix", "[QQMatrix]") {
    auto m = return_matrix();
    qqStdOut() << m;
}

TEST_CASE( "Matrix could be changed using row", "[QQMatrix]" ) {
    qqStdOut()<<"Matrix could be changed using row\n";

    QQMatrix<double> m(2,2, {
        1., 2.,
        3., 4.});
    QQMatrix<double> mc = m;
    REQUIRE(mc==m);
    mc.set(0, 0, 10.);
    REQUIRE(mc(0,0)==10.);
    REQUIRE(m(0,0)==1.);
    qqStdOut() << m;
    m.setrow(0, qqmat<2>({5., 6.}));
    qqStdOut() << "changed =\n" << m;
    REQUIRE( m == qqmat<2>({
        5., 6.,
        3., 4.
    }));
    qqStdOut() << "row(1) = " << m.row(1);
    REQUIRE( m.row(1) == qqmat<2>({3., 4.}));
}
