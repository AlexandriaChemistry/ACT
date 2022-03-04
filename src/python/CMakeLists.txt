#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
# Installation path for ACT python code. Depends from CMAKE_INSTALL_PREFIX.
set(ACT_PYTHON_LIB_PATH ${CMAKE_INSTALL_PREFIX}/lib/site-packages/act)

# Configure 'sys.path' in the executable scripts, so they may find the library.
configure_file(gauss2molprop.cmakein gauss2molprop @ONLY)
configure_file(generate_ff.cmakein generate_ff @ONLY)
configure_file(generate_mp.cmakein generate_mp @ONLY)
configure_file(mol_csv_api.py.cmakein mol_csv_api.py @ONLY)
configure_file(molselect.cmakein molselect @ONLY)

# Library scripts
install(FILES atomic_heat_of_formation.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES dofit.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES elements.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES gaff_to_alexandria.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES get_csv_rows.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES get_mol_dict.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES jacobi.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES molprops.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES molutils.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES read_gaussian_log.py DESTINATION ${ACT_PYTHON_LIB_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mol_csv_api.py DESTINATION ${ACT_PYTHON_LIB_PATH})

# Executable scripts
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gauss2molprop DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/generate_ff DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/generate_mp DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/molselect DESTINATION bin)
install(PROGRAMS install_act DESTINATION bin)
install(PROGRAMS viewxvg DESTINATION bin)
install(PROGRAMS view_fitness DESTINATION bin)

