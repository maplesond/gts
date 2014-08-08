//  ********************************************************************
//  This file is part of GTS (Good Transcript Selector).
//
//  GTS is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portculis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portculis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************


#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE GTS
#endif
#include <boost/test/unit_test.hpp>

#include <boost/filesystem.hpp>

#include <genbank.hpp>

using std::cout;
using std::endl;

using gts::gb::Genbank;

BOOST_AUTO_TEST_SUITE(genbank)

BOOST_AUTO_TEST_CASE(load) {
    
    std::vector<shared_ptr<Genbank> > genbank;
    Genbank::load("resources/test.gb", genbank);
    BOOST_CHECK(genbank.size() == 2);

    Genbank::save("resources/test_make.gb", genbank);
}

BOOST_AUTO_TEST_SUITE_END()
