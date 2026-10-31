#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
# Installation path for ACT python code. Depends from CMAKE_INSTALL_PREFIX.
set(ACT_PYTHON_LIB_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_PYTHON_LIBDIR}/act)

# Configure 'sys.path' in the executable scripts, so they may find the library.
configure_file(gauss2molprop.cmakein gauss2molprop @ONLY)
configure_file(generate_mp.cmakein generate_mp @ONLY)
configure_file(mol_csv_api.py.cmakein mol_csv_api.py @ONLY)
configure_file(molselect.cmakein molselect @ONLY)
configure_file(reshuffle_selection.cmakein reshuffle_selection @ONLY)

# Library scripts
install(FILES act.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES act_gct.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES act_openmm.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES actutils.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES atomic_heat_of_formation.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES dofit.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES elements.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES get_csv_rows.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES get_mol_dict.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES jacobi.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES molprops.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES molutils.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES ncia_reader_from_xyz.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES read_gaussian_log.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mol_csv_api.py DESTINATION ${ACT_PYTHON_LIB_PATH})

# Executable scripts
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gauss2molprop DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/generate_mp DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/molselect DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/reshuffle_selection DESTINATION bin)
install(PROGRAMS coords2molprop DESTINATION bin)
install(PROGRAMS dimer_scan DESTINATION bin)
install(PROGRAMS plot_convergence DESTINATION bin)
install(PROGRAMS view_fitness DESTINATION bin)
install(PROGRAMS donchev2molprop DESTINATION bin)
install(PROGRAMS ncia2molprop DESTINATION bin)
install(PROGRAMS actgui DESTINATION bin)

# Tests of the python scripts, they are skipped if OpenMM or RDKit cannot be imported
if (BUILD_TESTING AND GMX_BUILD_UNITTESTS)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if (Python3_Interpreter_FOUND)
        add_test(NAME ActOpenMMTest
                 COMMAND ${Python3_EXECUTABLE} -m unittest -v test_act_openmm
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
        set_tests_properties(ActOpenMMTest PROPERTIES
                             ENVIRONMENT "ACT_ALEXANDRIA=$<TARGET_FILE:alexandria>")
        # Skipped if RDKit cannot be imported
        add_test(NAME Coords2MolpropTest
                 COMMAND ${Python3_EXECUTABLE} -m unittest -v test_coords2molprop
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    endif()
endif()
