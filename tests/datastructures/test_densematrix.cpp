#include <catch2/catch.hpp>

#include <densematrix.h>

TEST_CASE("dense matrices", "[matrix]") {
        SECTION("matrix initialization"){
                densematrix m = densematrix();
                REQUIRE(m.countrows() == 0);
                REQUIRE(m.countcolumns() == 0);
        }
}
