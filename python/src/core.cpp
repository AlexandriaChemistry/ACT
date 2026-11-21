/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * The Python module act._core: the classes of ACT in Python.
 *
 * The package act runs the alexandria program for whole calculations. This
 * module gives Python the classes that the program is made of, such that a
 * script can read a force field or a molecule, look at it, change it and
 * compute with it, without a program in between. It is built with
 * -DACT_PYTHON=ON and needs pybind11.
 *
 * ACT is a library here: an error must not end Python. When the module is
 * loaded it lets fatal errors throw, see gmx_set_fatal_error_throws(), and
 * every exception of ACT becomes one of the package act: ActFileError,
 * ActInputError, ActInternalError or their base class ActError.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

// Before the headers of ACT, see there
#include "act_pybind.h"

#include <exception>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

namespace py = pybind11;

namespace alexandria
{

namespace
{

/*! \brief Raise an exception of the Python package act
 * \param[in] name    The name of the class in act.errors
 * \param[in] message The message
 */
void raise(const char *name, const char *message)
{
    // The classes are those of the package, such that one except clause
    // catches the errors of the program and of the classes alike
    const py::object cls = py::module_::import("act.errors").attr(name);
    PyErr_SetString(cls.ptr(), message);
}

/*! \brief Turn the exceptions of ACT into those of the Python package act
 *
 * By the error code of the exception, not by its class: this module and
 * the libraries of ACT each have their own copy of the type information of
 * the classes that are defined in a header only, such as
 * gmx::UserInputError, and a catch clause for such a class does not catch
 * an exception that another library threw.
 * \param[in] pointer The exception
 */
void translateException(std::exception_ptr pointer)
{
    try
    {
        if (pointer)
        {
            std::rethrow_exception(pointer);
        }
    }
    catch (const gmx::GromacsException &ex)
    {
        switch (ex.errorCode())
        {
            case gmx::eeFileNotFound:
            case gmx::eeFileIO:
                raise("ActFileError", ex.what());
                break;
            case gmx::eeInvalidInput:
            case gmx::eeInconsistentInput:
            case gmx::eeInvalidValue:
            case gmx::eeRange:
                raise("ActInputError", ex.what());
                break;
            case gmx::eeInternalError:
            case gmx::eeAPIError:
            case gmx::eeInvalidCall:
            case gmx::eeNotImplemented:
                raise("ActInternalError", ex.what());
                break;
            default:
                // A fatal error without a kind, a tolerance that is not met, and so on
                raise("ActError", ex.what());
                break;
        }
    }
}

/*! \brief Read a force field file
 * \param[in] fileName The file
 * \return The force field
 * \throws gmx::FileIOError if the file cannot be read, gmx::InvalidInputError
 * if what is in it is not accepted
 */
std::unique_ptr<ForceField> readForceFieldFile(const std::filesystem::path &fileName)
{
    auto pd = std::make_unique<ForceField>();
    readForceField(fileName.string(), pd.get());
    // As it was given: reading keeps the file it found, which for a name
    // without directory may be one of the library of ACT
    return pd;
}

//! \return The names of the particle types of a force field, in its order
std::vector<std::string> particleTypeNames(const ForceField &pd)
{
    std::vector<std::string> names;
    for (const auto &ptype : pd.particleTypesConst())
    {
        names.push_back(ptype.first.id());
    }
    return names;
}

//! \return The names of the interactions that a force field has, in its order
std::vector<std::string> interactionNames(const ForceField &pd)
{
    std::vector<std::string> names;
    for (const auto &force : pd.forcesConst())
    {
        names.push_back(interactionTypeToString(force.first));
    }
    return names;
}

//! Tell Python about the class ForceField
void defineForceField(py::module_ &module)
{
    py::class_<ForceField>(module, "ForceField",
                           "A force field: its particle types and the parameters of its interactions.\n\n"
                           "Read one from a file with ForceField.read(file).")
        .def_static("read", &readForceFieldFile, py::arg("file"),
                    "Read a force field from an XML file.\n\n"
                    "A name without a directory is also looked for in the library of ACT.\n"
                    "Raises ActFileError if the file cannot be read, and ActInputError if\n"
                    "what is in it is not a consistent force field.")
        .def("write",
             [](const ForceField &pd, const std::filesystem::path &fileName, bool compress)
             {
                 writeForceField(fileName.string(), &pd, compress);
             },
             py::arg("file"), py::arg("compress") = false,
             "Write the force field to an XML file, compressed if asked.\n\n"
             "Raises ActFileError if the file cannot be written.")
        .def_property_readonly("filename", [](const ForceField &pd) { return pd.filename(); },
                               "The file the force field was read from.")
        .def_property_readonly("polarizable", &ForceField::polarizable,
                               "Whether the force field has particles that are polarizable.")
        .def_property_readonly("particle_types", &particleTypeNames,
                               "The names of the particle types: atoms, shells and virtual sites.")
        .def_property_readonly("interactions", &interactionNames,
                               "The names of the interactions that the force field has parameters for.")
        .def("__repr__",
             [](const ForceField &pd)
             {
                 return "<act.ForceField " + pd.filename() + ": " +
                        std::to_string(pd.particleTypesConst().size()) + " particle types, " +
                        std::to_string(pd.forcesConst().size()) + " interactions>";
             });
}

} // namespace

} // namespace alexandria

PYBIND11_MODULE(_core, module)
{
    module.doc() = "The classes of ACT in Python. Use them through the package act.";
    // A library must leave its caller alive: errors are exceptions from here on
    gmx_set_fatal_error_throws(true);
    py::register_exception_translator(&alexandria::translateException);
    alexandria::defineForceField(module);
}
