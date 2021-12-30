#
#
#  This source file is part of the Alexandria program.
#
#  Copyright (C) 2014-2021
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

set(ALEXANDRIA_SOURCES
           acmfitnesscomputer.cpp
           acmindividual.cpp
           acminitializer.cpp
           alex_modules.cpp
      	   molpropobservable.cpp
           experiment.cpp
      	   molprop.cpp
           allbondeds.cpp
           analyze.cpp
           atype_mapping.cpp
           babel_io.cpp
           bastat.cpp
           bayes.cpp
           categories.cpp
           chargemodel.cpp
           communication.cpp
           composition.cpp
           confighandler.cpp
           devcomputer.cpp
           dissociation_energy.cpp
           fill_inputrec.cpp
           forcefieldparameter.cpp
           forcefieldparameterlist.cpp
#           gen_table.cpp
           gentop.cpp
 #          gentop_vsite.cpp
           gmx_simple_comm.cpp
           identifier.cpp
           interactiontype.cpp
           latex_util.cpp
           mcmcmutator.cpp
           memory_check.cpp
           merge_mp.cpp
           merge_pd.cpp
           molgen.cpp
           molprop_check.cpp
           molprop_csv.cpp
           molprop_sqlite3.cpp
           molprop_tables.cpp
           molprop_test.cpp
           molprop_util.cpp
           molprop_xml.cpp
           molselect.cpp
           mutability.cpp
           mymol.cpp
           mymol_low.cpp
           npointcrossover.cpp
           openmm_xml.cpp
           particletype.cpp
           percentmutator.cpp
           phase.cpp
  #         plistwrapper.cpp
           poldata.cpp
           poldata_low.cpp
           poldata_tables.cpp
           poldata_edit.cpp
           poldata_xml.cpp
           qgen_acm.cpp
           qgen_cube.cpp
           qgen_resp.cpp
           qm2molprop.cpp
           qtype.cpp
           readpsi4.cpp
           regression.cpp
           sharedindividualinfo.cpp
           stringutil.cpp
           symmetrize_charges.cpp
           topology.cpp
           tune_eem.cpp
#           tune_fc.cpp
#           tune_fc_utils.cpp
           tuning_utility.cpp          
           units.cpp
           xml_util.cpp
           analyze.cpp
           molprop_tables.cpp
           vsite.cpp
	)
	
add_subdirectory(coulombintegrals)

add_library(alexandria_objlib OBJECT ${ALEXANDRIA_SOURCES})

add_executable(alexandria alexandria.cpp
         $<TARGET_OBJECTS:alexandria_objlib>)

target_link_libraries(alexandria libactgromacs actga ${GMX_EXTRA_LIBRARIES} ${XML_LIBRARIES}
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

