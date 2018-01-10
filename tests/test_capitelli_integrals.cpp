/*
 * Copyright 2017 von Karman Institute for Fluid Dynamics (VKI)
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
#include <catch/catch.hpp>
#include <eigen3/Eigen/Dense>

using namespace Mutation;
using namespace Mutation::Utilities;
using namespace Mutation::Utilities::IO;
using namespace Catch;
using namespace Eigen;




/**
 * Tests the XML classes
 */
TEST_CASE
(
    "Laricchuita integrals",
    "[transport]"
)
{
    TemporaryFile file;
    file << "<tag att=\"100\">\n"
         << "    <!-- comment -->\n"
         << "    <child1> text </child1>\n"
         << "    <child2> <child3 att=\"value\"/> </child2>\n"
         << "</tag>";
    file.close();

    XmlDocument doc(file.filename());
    CHECK(doc.file() == file.filename());

    XmlElement& root = doc.root();
    CHECK(root.tag() == "tag");
    CHECK(root.text().empty());
    CHECK(root.line() == 1);
    CHECK(root.document() == &doc);

    std::string att_str; int att_int;
    CHECK(root.hasAttribute("att"));
    CHECK(root.getAttribute("att", att_str) == "100");
    CHECK(root.getAttribute("att", att_int) == 100);
    CHECK_FALSE(root.hasAttribute("bla"));

    XmlElement::const_iterator it = root.begin();
    CHECK(it->tag() == "child1");
    CHECK(it->text() == " text ");
    CHECK(it->line() == 3);
    CHECK(it->begin() == it->end());

    it++;
    CHECK(it->tag() == "child2");
    CHECK(it->text().empty());
    CHECK(it->line() == 4);
    CHECK(it+1 == root.end());

    it = it->begin();
    CHECK(it->tag() == "child3");
    CHECK(it->text().empty());
    CHECK(it->line() == 4);
    CHECK(it->getAttribute("att", att_str) == "value");

    // Parsing errors
    CHECK_THROWS_AS(XmlElement("a <tag/>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag"), Error);
    CHECK_THROWS_AS(XmlElement("<tag>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag att=\" />"), Error);
    CHECK_THROWS_AS(XmlElement("<tag att=\"\" >"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> text <tag>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> text </ta>"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> text </tag"), Error);
    CHECK_THROWS_AS(XmlElement("<tag> <child> </tag>"), Error);

}

