#include <catch.hpp>

#include <Constants.h>
#include <FourVector.h>
#include <ThreeVector.h>
#include <logging.h>

#include <cmath>

TEST_CASE( "Vector" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    SECTION( "Initialization" ) {
        yap::Vector<double, 3> v = {};
        yap::Vector<double, 3> v0({0, 0, 0});
        REQUIRE( v == v0);
    }

    SECTION( "ThreeVector" ) {

        auto zero = yap::ThreeVector<double>({0, 0, 0});
        auto v1 = yap::ThreeVector<double>({1, 2, 3});
        auto v2 = yap::ThreeVector<double>({4, 5, 6});

        SECTION( "addition-assignment" ) {
            v1 += v2;
            REQUIRE( v1 == yap::ThreeVector<double>({5, 7, 9}) );
        }

        SECTION( "subtraction-assignment" ) {
            v1 -= v2;
            REQUIRE( v1 == yap::ThreeVector<double>({ -3, -3, -3}) );
        }

        SECTION( "multiplication-assignment" ) {
            v1 *= 3.;
            REQUIRE( v1 == yap::ThreeVector<double>({3, 6, 9}) );
        }

        SECTION( "arithmetic operations" ) {

            // +
            REQUIRE( v1 + v2 == yap::ThreeVector<double>({5, 7, 9}) );

            // -
            REQUIRE( v1 - v2 == yap::ThreeVector<double>({ -3, -3, -3}) );

            // inner product
            REQUIRE( v1 * v2 == 32 );

            // norm
            REQUIRE( norm(v1) == 14 );
            REQUIRE( norm(v2) == 77 );

            // abs
            REQUIRE( abs(v1) == sqrt(14) );
            REQUIRE( abs(v2) == sqrt(77) );

            // cross product
            auto x = yap::ThreeVector<double>({1, 0, 0});
            auto y = yap::ThreeVector<double>({0, 1, 0});
            auto z = yap::ThreeVector<double>({0, 0, 1});
            REQUIRE( cross(x, y) == z );

            // (unary) minus
            REQUIRE( -v1 == -1. * v1 );
            REQUIRE( -v1 + v1 == zero );


        }

		SECTION( "commutativity" ) {
			auto A = yap::Vector<double, 4>({1,2,3,4});
			double a = 2;

			REQUIRE( a * A == A * a );

			auto B = yap::FourVector<double>({1,2,3,4});
			double b = 2;

			REQUIRE( b * B == B * b );
			std::cout << to_string( static_cast<yap::FourVector<double>>(b*B) ) << std::endl;

			auto C = yap::ThreeVector<double>({1,2,3});
			double c = 2;

			REQUIRE( c * C == C * c );
			std::cout << to_string( static_cast<yap::ThreeVector<double>>(c*C) ) << std::endl;

		}

        SECTION( "angles" ) {

            auto a = yap::ThreeVector<double>({2, 0, 0});
            auto b = yap::ThreeVector<double>({0, 2, 0});
            auto c = yap::ThreeVector<double>({0, 0, 2});
            auto z = yap::ThreeVector<double>({0, 0, 0});

            REQUIRE(angle(a, a) == 0.);
            REQUIRE(angle(a, -a) == yap::pi<double>());
            REQUIRE(angle(a, b) == 0.5 * yap::pi<double>());
            REQUIRE(angle(a, c) == 0.5 * yap::pi<double>());
            REQUIRE(std::isnan(angle(a, z)));

            REQUIRE( yap::theta(c, yap::ThreeAxes) == Approx(0.) );
            REQUIRE( yap::theta(a, yap::ThreeAxes) == Approx(yap::pi<double>() / 2.) );
            REQUIRE( yap::theta(b, yap::ThreeAxes) == Approx(yap::pi<double>() / 2.) );

            REQUIRE( yap::phi(a, yap::ThreeAxes) == Approx(0.) );
            REQUIRE( yap::phi(b, yap::ThreeAxes) == Approx(yap::pi<double>() / 2.) );
        }

        SECTION( "cross product" ) {
            REQUIRE( cross(yap::ThreeAxis_X, yap::ThreeAxis_Y) == yap::ThreeAxis_Z);
            REQUIRE( cross(yap::ThreeAxis_Y, yap::ThreeAxis_X) == -yap::ThreeAxis_Z);

            REQUIRE( cross(yap::ThreeAxis_Y, yap::ThreeAxis_Z) == yap::ThreeAxis_X);
            REQUIRE( cross(yap::ThreeAxis_Z, yap::ThreeAxis_Y) == -yap::ThreeAxis_X);

            REQUIRE( cross(yap::ThreeAxis_Z, yap::ThreeAxis_X) == yap::ThreeAxis_Y);
            REQUIRE( cross(yap::ThreeAxis_X, yap::ThreeAxis_Z) == -yap::ThreeAxis_Y);
        }

        SECTION( "constants" ) {

            // axes:
            REQUIRE( cross(yap::ThreeAxis_X, yap::ThreeAxis_Y) == yap::ThreeAxis_Z );
            REQUIRE( isRightHanded(yap::ThreeAxes) );
        }

    }

    SECTION( "FourVector" ) {

        auto v1 = yap::FourVector<double>({4, 3, 2, 1});
        auto v2 = yap::FourVector<double>({8, 7, 6, 5});

        SECTION( "addition-assignment" ) {
            v1 += v2;
            REQUIRE( v1 == yap::FourVector<double>({12, 10, 8, 6}) );
        }

        SECTION( "subtraction-assignment" ) {
            v1 -= v2;
            REQUIRE( v1 == yap::FourVector<double>({ -4, -4, -4, -4}) );
        }

        SECTION( "multiplication-assignment" ) {
            v1 *= 2.;
            REQUIRE( v1 == yap::FourVector<double>({8, 6, 4, 2}) );
        }

        SECTION( "arithmetic operations" ) {

            // +
            REQUIRE( v1 + v2 == yap::FourVector<double>({12, 10, 8, 6}) );

            // -
            REQUIRE( v1 - v2 == yap::FourVector<double>({ -4, -4, -4, -4}) );

            // * (four-vector inner product)
            REQUIRE( v1 * v2 == -6 );

            // norm
            REQUIRE( norm(v1) == 2 );
            REQUIRE( norm(v2) == -46 );

            // abs
            REQUIRE( abs(v1) == sqrt(2) );
            REQUIRE( abs(yap::FourVector<double>({1, 0, 0, 0})) == 1 );
            REQUIRE( abs(yap::FourVector<double>({1, 1, 0, 0})) == 0 );

        }

        // unit
        //REQUIRE( unit(yap::FourVector<double>({1,0,0,0})) == yap::FourVector<double>({1,0,0,0}) );

        // *

        // -
        SECTION("minus") {
            REQUIRE( -v1 == yap::FourVector<double>({4, -3, -2, -1}) );

            std::vector<yap::FourVector<double>> VV = {v1, v2};
            std::vector<yap::FourVector<double>> mVV = { -v1, -v2};
            REQUIRE( -VV == mVV );
        }

        // cross

        // angle

    }

}
