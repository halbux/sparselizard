#include <catch2/catch.hpp>

#include <intdensematrix.h>

TEST_CASE("integer dense matrices", "[int matrix]") {
        SECTION("matrix initialization"){
                intdensematrix m = intdensematrix();
                REQUIRE(m.countrows() == 0);
                REQUIRE(m.countcolumns() == 0);
        }
}
