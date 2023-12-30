
file(GLOB STLBFGS_SOURCES *.cpp *.h)
target_sources(actLib PRIVATE ${STLBFGS_SOURCES})
	
if (BUILD_TESTING)
#    add_subdirectory(tests)
endif()

