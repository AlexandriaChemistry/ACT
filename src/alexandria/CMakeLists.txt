#
#
#  This source file is part of the Alexandria program.
#
#  Copyright (C) 2014-2022
#
#  Developers:
#             Mohammad Mehdi Ghahremanpour,
#             Paul J. van Maaren,
#             David van der Spoel (Project leader)
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.
#

target_sources(actLib
    PRIVATE
           acm_ga.cpp
           acmfitnesscomputer.cpp
           acmindividual.cpp
           acminitializer.cpp
           acthelper.cpp
           actmiddleman.cpp
           alex_modules.cpp
           allbondeds.cpp
           analyze.cpp
           analyze.cpp
           atype_mapping.cpp
           babel_io.cpp
           bastat.cpp
           bayes.cpp
           chargemodel.cpp
           composition.cpp
           confighandler.cpp
           devcomputer.cpp
           dissociation_energy.cpp
           fill_inputrec.cpp
           gentop.cpp
           gromacs_top.cpp
           mcmcmutator.cpp
           merge_mp.cpp
           merge_pd.cpp
           molgen.cpp
           molselect.cpp
           mymol.cpp
           mymol_low.cpp
           npointcrossover.cpp
           openmm_xml.cpp
           percentmutator.cpp
           qgen_acm.cpp
           qgen_cube.cpp
           qgen_resp.cpp
           qm2molprop.cpp
           qtype.cpp
           readpsi4.cpp
           staticindividualinfo.cpp
           symmetrize_charges.cpp
           topology.cpp
           tune_ff.cpp
           tuning_utility.cpp          
           vsite.cpp
	)
	
add_executable(alexandria alexandria.cpp ${actLib})

target_link_libraries(alexandria actLib libactgromacs 
                      ${GMX_EXTRA_LIBRARIES} ${XML_LIBRARIES}
                      ${GMX_EXE_LINKER_FLAGS})
set(BINARY_NAME "alexandria${GMX_BINARY_SUFFIX}")
set_target_properties(alexandria PROPERTIES 
                      OUTPUT_NAME "${BINARY_NAME}")

install(TARGETS alexandria
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libalexandria.pc.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libalexandria.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
        RENAME "libalexandria${GMX_LIBS_SUFFIX}.pc"
        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
