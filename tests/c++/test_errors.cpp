/*
 * Copyright 2017-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "mutation++.h"
#include <catch.hpp>

using namespace Mutation;
using namespace Catch;

/**
 * Runs basic tests on the initialization and use of the Mutation++ exception
 * classes defined in Errors.h.
 */
TEST_CASE("Test exception classes", "[exceptions]"
)
{
    SECTION("Error (base class)") {
        Error error = Error("type")("info1", "value") << "message";
        error("info2", 37);
        error.addExtraInfo("info3", 'a');
        CHECK(error.getExtraInfo("info1") == "value");
        CHECK(error.getExtraInfo("info2") == "37");
        CHECK(error.getExtraInfo("info3") == "a");
        CHECK(error.getExtraInfo("bla") == "");

        CHECK(
            std::string(error.what()) ==
            "\nM++ error: type.\ninfo1: value\ninfo2: 37\ninfo3: a\nmessage\n"
        );

        CHECK_THROWS_AS(throw Error("type"), Error);
    }

    SECTION("FileNotFoundError") {
        FileNotFoundError error = FileNotFoundError("name") << "message";
        CHECK(error.getExtraInfo("file") == "name");
        CHECK(error.file() == "name");

        CHECK(
            std::string(error.what()) ==
            (Error("file not found")("file", "name") << "message").what()
        );

        CHECK_THROWS_AS(throw FileNotFoundError("") << "", Error);
        CHECK_THROWS_AS(throw FileNotFoundError("") << "", FileNotFoundError);
    }

    SECTION("FileParseError") {
        FileParseError error = FileParseError("name", 18) << "message";
        CHECK(error.getExtraInfo("file") == "name");
        CHECK(error.file() == "name");
        CHECK(error.getExtraInfo("line") == "18");
        CHECK(error.line() == "18");

        CHECK(
            std::string(error.what()) ==
            (Error("error parsing file")("file", "name")("line", 18) << "message").what()
        );

        CHECK_THROWS_AS(throw FileParseError("", 0) << "", Error);
        CHECK_THROWS_AS(throw FileParseError("", 0) << "", FileParseError);
    }

    SECTION("InvalidInputError") {
        InvalidInputError error = InvalidInputError("input", 37) << "message";
        CHECK(error.getExtraInfo("input") == "input = 37");
        CHECK(error.inputName() == "input");
        CHECK(error.inputValue() == "37");

        CHECK(
            std::string(error.what()) ==
            (Error("invalid input")("input", "input = 37") << "message").what()
        );

        CHECK_THROWS_AS(throw InvalidInputError("", 0) << "", Error);
        CHECK_THROWS_AS(throw InvalidInputError("", 0) << "", InvalidInputError);
    }

    SECTION("LogicError") {
        std::stringstream ss;
        LogicError error = LogicError() << "message"; ss << __LINE__;
        CHECK(error.getExtraInfo("file") == __FILE__);
        CHECK(error.file() == __FILE__);
        CHECK(error.getExtraInfo("line") == ss.str());
        CHECK(error.line() == ss.str());

        CHECK(
            std::string(error.what()) ==
            (Error("logic error")("file", __FILE__)("line", ss.str()) << "message").what()
        );

        CHECK_THROWS_AS(throw LogicError() << "", Error);
        CHECK_THROWS_AS(throw LogicError() << "", LogicError);
    }

    SECTION("NotImplementedError") {
        std::stringstream ss;
        NotImplementedError error = NotImplementedError("function") << "message"; ss << __LINE__;
        CHECK(error.getExtraInfo("file") == __FILE__);
        CHECK(error.file() == __FILE__);
        CHECK(error.getExtraInfo("line") == ss.str());
        CHECK(error.line() == ss.str());
        CHECK(error.getExtraInfo("function") == "function");
        CHECK(error.function() == "function");

        CHECK(
            std::string(error.what()) ==
            (Error("function not implemented")("function", "function")
                  ("file", __FILE__)("line", ss.str()) << "message").what()
        );

        CHECK_THROWS_AS(throw NotImplementedError("") << "", Error);
        CHECK_THROWS_AS(throw NotImplementedError("") << "", NotImplementedError);
    }
}
