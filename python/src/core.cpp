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

#include <cmath>
#include <exception>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_json.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_validation.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/particletype.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

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

/*! \brief A parameter of a force field as Python sees it
 * What the parameter was when it was asked for. To change it, see
 * setParameter: a copy cannot outlive the force field or go out of step
 * with it.
 */
struct Parameter
{
    //! The interaction it belongs to, empty for a parameter of a particle type
    std::string interaction;
    //! The identifier of the types it is for, or the name of the particle type
    std::string identifier;
    //! The name of the parameter
    std::string name;
    //! The value, its bounds and its uncertainty, in the unit
    double      value = 0, minimum = 0, maximum = 0, uncertainty = 0;
    //! The unit
    std::string unit;
    //! Whether training changes it: Bounded, Fixed, ACM or Dependent
    std::string mutability;
    //! The number of data points it was trained on
    int         ntrain = 0;
    //! Whether it may not be negative
    bool        nonnegative = false;
};

//! \return A parameter for Python
Parameter makeParameter(const std::string &interaction, const std::string &identifier,
                        const std::string &name, const ForceFieldParameter &param)
{
    Parameter made;
    made.interaction = interaction;
    made.identifier  = identifier;
    made.name        = name;
    made.value       = param.value();
    made.minimum     = param.minimum();
    made.maximum     = param.maximum();
    made.uncertainty = param.uncertainty();
    made.unit        = param.unit();
    made.mutability  = mutabilityName(param.mutability());
    made.ntrain      = static_cast<int>(param.ntrain());
    made.nonnegative = param.nonNegative();
    return made;
}

/*! \brief Look up an interaction by its name in the force field file
 * \param[in] name The name, such as VANDERWAALS
 * \return The interaction
 * \throws gmx::InvalidInputError if ACT has no interaction of that name
 */
InteractionType interactionOf(const std::string &name)
{
    InteractionType itype;
    if (!stringToInteractionType(name, &itype))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("ACT has no interaction '%s'. The interactions of a force field are in its property interactions", name.c_str())));
    }
    return itype;
}

/*! \brief The parameters of an interaction
 * \param[in] pd          The force field
 * \param[in] interaction The name of the interaction
 * \return The parameters of all its entries, none if the force field does
 * not have the interaction
 * \throws gmx::InvalidInputError if ACT has no interaction of that name
 */
std::vector<Parameter> parametersOf(const ForceField &pd, const std::string &interaction)
{
    std::vector<Parameter> all;
    const auto             itype = interactionOf(interaction);
    if (!pd.interactionPresent(itype))
    {
        return all;
    }
    for (const auto &entry : pd.findForcesConst(itype).parametersConst())
    {
        for (const auto &param : entry.second)
        {
            all.push_back(makeParameter(interaction, entry.first.id(), param.first, param.second));
        }
    }
    return all;
}

/*! \brief Find a parameter of an interaction to look at or to change
 * \param[in] pd          The force field
 * \param[in] interaction The name of the interaction
 * \param[in] identifier  The identifier of the entry, as in the force field file
 * \param[in] name        The name of the parameter
 * \return The parameter
 * \throws gmx::InvalidInputError if there is no such parameter
 */
ForceFieldParameter *findParameter(ForceField *pd, const std::string &interaction,
                                   const std::string &identifier, const std::string &name)
{
    const auto itype = interactionOf(interaction);
    if (!pd->interactionPresent(itype))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The force field has no interaction %s", interaction.c_str())));
    }
    auto *forces = pd->findForces(itype);
    for (const auto &entry : forces->parametersConst())
    {
        if (entry.first.id() == identifier)
        {
            auto *params = forces->findParameters(entry.first);
            auto  found  = params->find(name);
            if (found == params->end())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("%s of %s has no parameter '%s'",
                                                                   identifier.c_str(), interaction.c_str(), name.c_str())));
            }
            return &found->second;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("%s has no entry '%s'", interaction.c_str(), identifier.c_str())));
}

/*! \brief Find a particle type
 * \param[in] pd   The force field
 * \param[in] name The name of the particle type
 * \return The particle type
 * \throws gmx::InvalidInputError if there is no such particle type
 */
ParticleType *findParticle(ForceField *pd, const std::string &name)
{
    for (auto &ptype : *pd->particleTypes())
    {
        if (ptype.first.id() == name)
        {
            return &ptype.second;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("The force field has no particle type '%s'", name.c_str())));
}

//! A particle type of a force field as Python sees it, see Parameter
struct Particle
{
    //! The name
    std::string                        name;
    //! Atom, Shell or VSite
    std::string                        kind;
    //! What the force field says about it
    std::string                        description;
    //! The options, such as the element and the types it has for each interaction
    std::map<std::string, std::string> options;
    //! Its own parameters: the charge, the mass and so on
    std::vector<Parameter>             parameters;
};

//! \return A particle type for Python
Particle particleOf(ForceField &pd, const std::string &name)
{
    const auto *ptype = findParticle(&pd, name);
    Particle    made;
    made.name        = name;
    made.kind        = actParticleToString(ptype->apType());
    made.description = ptype->description();
    made.options     = ptype->optionsConst();
    for (const auto &param : ptype->parametersConst())
    {
        made.parameters.push_back(makeParameter("", name, param.first, param.second));
    }
    return made;
}

/*! \brief Change a parameter, after a check that what comes out makes sense
 *
 * Nothing is changed if the change is refused. A parameter that training
 * does not change (Fixed) is changed all the same: that is what the caller
 * asks for.
 * \param[in] param      The parameter
 * \param[in] what       What it is, for messages
 * \param[in] value      The new value, if any
 * \param[in] minimum    The new minimum, if any
 * \param[in] maximum    The new maximum, if any
 * \param[in] mutability The new mutability, if any
 * \throws gmx::InvalidInputError if nothing is to be changed, a number is
 * not finite, the mutability is unknown, the minimum would be above the
 * maximum, the value outside the bounds, or a parameter negative that
 * cannot be
 */
void changeParameter(ForceFieldParameter *param, const std::string &what,
                     std::optional<double> value, std::optional<double> minimum,
                     std::optional<double> maximum, const std::optional<std::string> &mutability)
{
    if (!value && !minimum && !maximum && !mutability)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Nothing to change for %s: give value, minimum, maximum or mutability", what.c_str())));
    }
    for (const auto &number : { value, minimum, maximum })
    {
        if (number && !std::isfinite(*number))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("A number for %s is not finite", what.c_str())));
        }
    }
    Mutability mut = param->mutability();
    if (mutability && !nameToMutability(*mutability, &mut))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown mutability '%s' for %s: use Bounded or Fixed", mutability->c_str(), what.c_str())));
    }
    const double newValue = value ? *value : param->value();
    const double newMin   = minimum ? *minimum : param->minimum();
    const double newMax   = maximum ? *maximum : param->maximum();
    if (newMin > newMax)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The minimum %g of %s would be above its maximum %g",
                                                           newMin, what.c_str(), newMax)));
    }
    if (newValue < newMin || newValue > newMax)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value %g of %s would be outside its bounds %g to %g: give bounds that include it",
                                                           newValue, what.c_str(), newMin, newMax)));
    }
    if (param->nonNegative() && newMin < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("%s cannot be negative, so its minimum cannot be %g", what.c_str(), newMin)));
    }
    // The bounds wide enough for the old and the new value on the way
    param->setMaximum(std::max(newMax, param->maximum()));
    (void) param->setMinimum(std::min(newMin, param->minimum()));
    param->forceSetValue(newValue);
    (void) param->setMinimum(newMin);
    param->setMaximum(newMax);
    param->setMutability(mut);
}

//! \return The problems with the parameters of a force field, for Python
py::list validationOf(const ForceField &pd)
{
    py::list issues;
    for (const auto &issue : validateForceField(pd))
    {
        py::dict one;
        one["kind"]       = forceFieldIssueKindName(issue.kind);
        one["where"]      = issue.where;
        one["identifier"] = issue.identifier;
        one["parameter"]  = issue.parameter;
        one["text"]       = issue.text;
        one["repaired"]   = issue.repaired;
        issues.append(one);
    }
    return issues;
}

//! Tell Python about the parameters and the particle types of a force field
void defineParameters(py::module_ &module)
{
    py::class_<Parameter>(module, "Parameter",
                          "A parameter of a force field, as it was when it was asked for.\n\n"
                          "Change it with ForceField.set_parameter or set_particle_parameter.")
        .def_readonly("interaction", &Parameter::interaction,
                      "The interaction it belongs to, empty for a parameter of a particle type.")
        .def_readonly("identifier", &Parameter::identifier,
                      "The types it is for, as in the force field file, or the particle type.")
        .def_readonly("name", &Parameter::name, "The name of the parameter, such as sigma.")
        .def_readonly("value", &Parameter::value, "The value, in the unit.")
        .def_readonly("minimum", &Parameter::minimum, "The smallest value that training gives it.")
        .def_readonly("maximum", &Parameter::maximum, "The largest value that training gives it.")
        .def_readonly("uncertainty", &Parameter::uncertainty, "The uncertainty of the value.")
        .def_readonly("unit", &Parameter::unit, "The unit.")
        .def_readonly("mutability", &Parameter::mutability,
                      "Whether training changes it: Bounded, Fixed, ACM or Dependent.")
        .def_readonly("ntrain", &Parameter::ntrain, "The number of data points it was trained on.")
        .def_readonly("nonnegative", &Parameter::nonnegative, "Whether it may not be negative.")
        .def("__repr__",
             [](const Parameter &param)
             {
                 return gmx::formatString("<act.Parameter %s%s%s %s = %g %s, %g to %g, %s>",
                                          param.interaction.c_str(), param.interaction.empty() ? "" : " ",
                                          param.identifier.c_str(), param.name.c_str(), param.value,
                                          param.unit.c_str(), param.minimum, param.maximum,
                                          param.mutability.c_str());
             });
    py::class_<Particle>(module, "ParticleType",
                         "A particle type of a force field, as it was when it was asked for.")
        .def_readonly("name", &Particle::name, "The name.")
        .def_readonly("kind", &Particle::kind, "Atom, Shell or VSite.")
        .def_readonly("description", &Particle::description, "What the force field says about it.")
        .def_readonly("options", &Particle::options,
                      "The options: the element, and the types it has for each interaction.")
        .def_readonly("parameters", &Particle::parameters,
                      "Its own parameters, such as the charge and the mass.")
        .def("__repr__",
             [](const Particle &particle)
             {
                 return gmx::formatString("<act.ParticleType %s (%s), %zu parameters>", particle.name.c_str(),
                                          particle.kind.c_str(), particle.parameters.size());
             });
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
             [](ForceField &pd, const std::filesystem::path &fileName, bool compress)
             {
                 // The time and the checksum of what is written now
                 pd.updateTimeStamp();
                 pd.updateCheckSum();
                 writeForceField(fileName.string(), &pd, compress);
             },
             py::arg("file"), py::arg("compress") = false,
             "Write the force field to an XML file, compressed if asked.\n\n"
             "Raises ActFileError if the file cannot be written.")
        .def("parameters", &parametersOf, py::arg("interaction"),
             "The parameters of an interaction, such as VANDERWAALS, for all its entries.\n\n"
             "An empty list if the force field does not have the interaction. Raises\n"
             "ActInputError if ACT has no interaction of that name.")
        .def("parameter",
             [](ForceField &pd, const std::string &interaction, const std::string &identifier,
                const std::string &name)
             {
                 return makeParameter(interaction, identifier, name,
                                      *findParameter(&pd, interaction, identifier, name));
             },
             py::arg("interaction"), py::arg("identifier"), py::arg("name"),
             "One parameter: of an interaction, for the types of the identifier as in the\n"
             "force field file (such as c3 or c3_b~h_b), by its name (such as sigma).\n\n"
             "Raises ActInputError if there is no such parameter.")
        .def("set_parameter",
             [](ForceField &pd, const std::string &interaction, const std::string &identifier,
                const std::string &name, std::optional<double> value, std::optional<double> minimum,
                std::optional<double> maximum, const std::optional<std::string> &mutability)
             {
                 changeParameter(findParameter(&pd, interaction, identifier, name),
                                 name + " of " + identifier + " (" + interaction + ")",
                                 value, minimum, maximum, mutability);
             },
             py::arg("interaction"), py::arg("identifier"), py::arg("name"), py::kw_only(),
             py::arg("value") = py::none(), py::arg("minimum") = py::none(),
             py::arg("maximum") = py::none(), py::arg("mutability") = py::none(),
             "Change the value, the bounds or the mutability (Bounded or Fixed) of a parameter.\n\n"
             "What is not given stays as it is. Raises ActInputError, and changes nothing, if\n"
             "the value would be outside the bounds, the minimum above the maximum, or a\n"
             "parameter negative that cannot be.")
        .def("particle_type", &particleOf, py::arg("name"),
             "A particle type with its kind, its options and its parameters.\n\n"
             "Raises ActInputError if the force field has no particle type of that name.")
        .def("set_particle_parameter",
             [](ForceField &pd, const std::string &particle, const std::string &name,
                std::optional<double> value, std::optional<double> minimum,
                std::optional<double> maximum, const std::optional<std::string> &mutability)
             {
                 auto *ptype = findParticle(&pd, particle);
                 if (!ptype->hasParameter(name))
                 {
                     GMX_THROW(gmx::InvalidInputError(gmx::formatString("The particle type %s has no parameter '%s'",
                                                                        particle.c_str(), name.c_str())));
                 }
                 changeParameter(ptype->parameter(name), name + " of particle type " + particle,
                                 value, minimum, maximum, mutability);
             },
             py::arg("particle"), py::arg("name"), py::kw_only(),
             py::arg("value") = py::none(), py::arg("minimum") = py::none(),
             py::arg("maximum") = py::none(), py::arg("mutability") = py::none(),
             "Change a parameter of a particle type, such as its charge, see set_parameter.")
        .def("validate", &validationOf,
             "The problems with the parameters, as edit_ff -validate lists them: an empty\n"
             "list if there are none, otherwise for each problem its kind, where it is and\n"
             "what is wrong.")
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
    alexandria::defineParameters(module);
    alexandria::defineForceField(module);
}
