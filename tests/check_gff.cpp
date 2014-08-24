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

#include <iostream>
using std::cout;
using std::endl;

#include <boost/filesystem.hpp>

#include <gff.hpp>
using gts::gff::GFFModel;

BOOST_AUTO_TEST_SUITE(gff)

BOOST_AUTO_TEST_CASE(gffLoad) {
    
    shared_ptr<GFFModel> geneModel = GFFModel::load("resources/test_tair10_head.gff");
    
    BOOST_CHECK(geneModel->getNbGenes() == 6);
    BOOST_CHECK(geneModel->getTotalNbTranscripts() == 8);
    BOOST_CHECK(1 == 1);
}


BOOST_AUTO_TEST_SUITE_END()
