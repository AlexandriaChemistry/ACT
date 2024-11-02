#
# This file is part of the Alexandria Chemistry Toolkit
#

set(testname "SobolTests")
set(exename "sobol-test")

gmx_add_unit_test(
    ${testname}
    ${exename}
    # files with code for tests
    sobolTest.cpp
    # pseudo-library for code for alexandria
    )
target_link_libraries(${exename} actLib libactgromacs ${GMX_EXTRA_LIBRARIES} ${XML_LIBRARIES} ${GMX_EXE_LINKER_FLAGS})
