#
# CMakeLists.txt for documentation generation using Sphinx
#
# This file builds both HTML website and PDF manual from reStructuredText sources
#

find_package(Sphinx REQUIRED)

# Set output directories
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
set(SPHINX_HTML_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
set(SPHINX_LATEX_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/latex")
set(SPHINX_PDF_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/pdf")
set(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Variables for Sphinx configuration
set(SPHINX_EXECUTABLE sphinx-build CACHE FILEPATH "Path to sphinx-build executable")

# Check if Sphinx build is possible
if(NOT SPHINX_EXECUTABLE)
    message(WARNING "Sphinx executable not found. Website and PDF manual cannot be built.")
    set(SPHINX_BUILD_POSSIBLE FALSE)
else()
    set(SPHINX_BUILD_POSSIBLE TRUE)
endif()

configure_file(conf-vars.py.cmakein ${SPHINX_CACHE_DIR}/conf-vars.py @ONLY)

if(SPHINX_BUILD_POSSIBLE)
    include(FindPythonModule)
    find_python_module(pygments ${_find_deps_options})
    message("Pygments_FOUND ${PYTHONMODULE_PYGMENTS}")
    if (PYTHONMODULE_PYGMENTS)
        set(Sphinx_pygments_FOUND 1)
    endif()
    
    # ============================================================================
    # HTML Website Generation
    # ============================================================================
    add_custom_command(
        OUTPUT "${SPHINX_HTML_OUTPUT_DIR}/index.html"
        DEPENDS ${SPHINX_SOURCE_DIR}/conf.py ${SPHINX_SOURCE_DIR}/index.rst
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b html
            -W --keep-going
            -w "${CMAKE_CURRENT_BINARY_DIR}/sphinx-html.log"
            -d "${SPHINX_CACHE_DIR}"
            "${SPHINX_SOURCE_DIR}"
            "${SPHINX_HTML_OUTPUT_DIR}"
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building HTML website documentation"
        VERBATIM
    )

    add_custom_target(website
        DEPENDS "${SPHINX_HTML_OUTPUT_DIR}/index.html"
        COMMENT "HTML website built in ${SPHINX_HTML_OUTPUT_DIR}"
    )

    # ============================================================================
    # LaTeX and PDF Manual Generation
    # ============================================================================
    add_custom_command(
        OUTPUT "${SPHINX_LATEX_OUTPUT_DIR}/act.tex"
        DEPENDS ${SPHINX_SOURCE_DIR}/conf.py ${SPHINX_SOURCE_DIR}/index.rst
        COMMAND
            ${SPHINX_EXECUTABLE}
            -q -b latex
            -W --keep-going
            -w "${CMAKE_CURRENT_BINARY_DIR}/sphinx-latex.log"
            -d "${SPHINX_CACHE_DIR}"
            "${SPHINX_SOURCE_DIR}"
            "${SPHINX_LATEX_OUTPUT_DIR}"
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building LaTeX intermediate files for PDF"
        VERBATIM
    )

    # Convert LaTeX to PDF using pdflatex
    find_program(PDFLATEX_EXECUTABLE pdflatex)
    
    if(PDFLATEX_EXECUTABLE)
        add_custom_command(
            OUTPUT "${SPHINX_PDF_OUTPUT_DIR}/act.pdf"
            DEPENDS "${SPHINX_LATEX_OUTPUT_DIR}/act.tex"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${SPHINX_PDF_OUTPUT_DIR}"
            COMMAND cd "${SPHINX_LATEX_OUTPUT_DIR}" && make
            COMMAND ${CMAKE_COMMAND} -E copy
                "${SPHINX_LATEX_OUTPUT_DIR}/act.pdf"
                "${SPHINX_PDF_OUTPUT_DIR}/act.pdf"
            WORKING_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Building PDF manual from LaTeX"
            VERBATIM
        )

        add_custom_target(manual
            DEPENDS "${SPHINX_PDF_OUTPUT_DIR}/act.pdf"
            COMMENT "PDF manual built in ${SPHINX_PDF_OUTPUT_DIR}"
        )
    else()
        message(WARNING "pdflatex not found. PDF manual cannot be built.")
        add_custom_target(manual
            COMMAND ${CMAKE_COMMAND} -E echo
            "Error: pdflatex not found. Cannot build PDF manual."
            VERBATIM
        )
    endif()

    # ============================================================================
    # Combined Documentation Target
    # ============================================================================
    if(PDFLATEX_EXECUTABLE)
        add_custom_target(docs
            DEPENDS website manual
            COMMENT "All documentation (website and manual) built successfully"
        )
    else()
        add_custom_target(docs
            DEPENDS website
            COMMENT "Website documentation built (PDF manual skipped - pdflatex not found)"
        )
    endif()

    # ============================================================================
    # Clean Target
    # ============================================================================
    add_custom_target(docs-clean
        COMMAND ${CMAKE_COMMAND} -E remove_directory "${SPHINX_HTML_OUTPUT_DIR}"
        COMMAND ${CMAKE_COMMAND} -E remove_directory "${SPHINX_LATEX_OUTPUT_DIR}"
        COMMAND ${CMAKE_COMMAND} -E remove_directory "${SPHINX_PDF_OUTPUT_DIR}"
        COMMAND ${CMAKE_COMMAND} -E remove_directory "${SPHINX_CACHE_DIR}"
        COMMENT "Cleaning up generated documentation"
        VERBATIM
    )

else()
    # Fallback targets when Sphinx is not available
    add_custom_target(website
        COMMAND ${CMAKE_COMMAND} -E echo
        "Cannot build website. Sphinx is not available."
        VERBATIM
    )

    add_custom_target(manual
        COMMAND ${CMAKE_COMMAND} -E echo
        "Cannot build PDF manual. Sphinx is not available."
        VERBATIM
    )

    add_custom_target(docs
        COMMAND ${CMAKE_COMMAND} -E echo
        "Cannot build documentation. Sphinx is not available."
        VERBATIM
    )

endif()
