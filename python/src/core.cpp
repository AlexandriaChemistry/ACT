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

#include <dlfcn.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/alexandria/train_utility.h"
#include "act/basics/atomization_energy.h"
#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/basics/mutability.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_json.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_validation.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/particletype.h"
#include "act/forces/forcecomputer.h"
#include "act/import/compound_reader.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_json.h"
#include "act/utility/jsontree.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/exceptions.h"
#include "buildinfo.h"
#include "gromacs/math/functions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

namespace py = pybind11;

namespace alexandria
{

namespace
{

/*! \brief Where ACT is installed, for the code that looks for its data files
 *
 * A program of ACT finds the data directory from where the program is. In
 * Python the program is the interpreter, which is somewhere else, so this
 * module tells where it is itself: in the package, a fixed number of
 * directories below the installation prefix. For a module in a build tree
 * it is the prefix that ACT was configured with.
 */
class ModuleContext : public gmx::IProgramContext
{
public:
    ModuleContext()
    {
        Dl_info info;
        // The file that this very function was loaded from
        if (dladdr(reinterpret_cast<const void *>(&anchor), &info) != 0 && info.dli_fname != nullptr)
        {
            std::error_code ignored;
            path_ = info.dli_fname;
            const auto installed = std::filesystem::weakly_canonical(std::filesystem::path(path_).parent_path() /
                                                                     ACT_PYTHON_TO_PREFIX, ignored);
            if (!ignored && std::filesystem::is_directory(installed / ACT_INSTALL_ACTDATADIR, ignored))
            {
                prefix_ = installed.string();
            }
        }
        std::error_code ignored;
        if (prefix_.empty() &&
            std::filesystem::is_directory(std::filesystem::path(CMAKE_INSTALL_PREFIX) / ACT_INSTALL_ACTDATADIR, ignored))
        {
            prefix_ = CMAKE_INSTALL_PREFIX;
        }
    }
    const char *programName() const override { return "act"; }
    const char *displayName() const override { return "The Python module of ACT"; }
    const char *fullBinaryPath() const override { return path_.c_str(); }
    gmx::InstallationPrefixInfo installationPrefix() const override
    {
        return gmx::InstallationPrefixInfo(prefix_.c_str(), false);
    }
    const char *commandLine() const override { return ""; }
    //! \return The data directory of the installation, empty if it was not found
    std::string dataDirectory() const
    {
        return prefix_.empty() ? std::string() : (std::filesystem::path(prefix_) / ACT_INSTALL_ACTDATADIR).string();
    }

private:
    //! A function of this module, to find the file it is in
    static void anchor() {}
    //! The file of this module and the installation prefix, empty if not found
    std::string path_, prefix_;
};

//! \return The context of this module, which lives as long as the process
const ModuleContext &moduleContext()
{
    static const ModuleContext context;
    return context;
}

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
    catch (const py::error_already_set &)
    {
        // An error of Python itself, which Python gets as it is
        throw;
    }
    catch (const py::builtin_exception &)
    {
        // A TypeError, ValueError and so on of this module
        throw;
    }
    catch (const std::bad_alloc &)
    {
        throw;
    }
    catch (const std::exception &ex)
    {
        // Of a library that ACT uses, such as the one that reads structure files
        raise("ActError", ex.what());
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

/*! \brief Give Python what a JSON tree has
 * The descriptions that the programs of ACT write in their json files are
 * what the classes tell Python, such that both say the same.
 * \param[in] tree The tree
 * \return A dict, a list, a number or a text
 */
py::object toPython(const JsonTree &tree)
{
    JsonTree root("root");
    root.addObject(tree);
    int        indent = 0;
    const auto text   = root.writeString(true, &indent);
    return py::module_::import("json").attr("loads")(text)["root"][py::str(tree.key())];
}

/*! \brief A message handler that keeps what it is told instead of printing it
 * A function of ACT that cannot do what it is asked tells its message
 * handler and returns. The messages are for the exception that Python gets.
 */
class QuietMessages
{
public:
    QuietMessages()
    {
        // Only what ends a program is printed, and that throws here
        handler_.setPrintLevel(ACTStatus::Fatal);
        handler_.record(&messages_, ACTStatus::Warning);
    }
    //! \return The message handler to give to a function of ACT
    MsgHandler *handler() { return &handler_; }
    /*! \brief Throw if a function reported an error
     * \param[in] what What was being done, for the message
     * \throws gmx::InvalidInputError with the errors that were reported
     */
    void check(const std::string &what) const
    {
        if (handler_.ok())
        {
            return;
        }
        std::string text;
        for (const auto &message : messages_)
        {
            if (message.level <= ACTStatus::Error && !message.text.empty())
            {
                text += (text.empty() ? "" : " ") + message.text;
            }
        }
        // Messages end with a new line of their own
        while (!text.empty() && (text.back() == '\n' || text.back() == ' '))
        {
            text.pop_back();
        }
        GMX_THROW(gmx::InvalidInputError(what + (text.empty() ? "" : ": " + text)));
    }

private:
    //! The message handler
    MsgHandler                      handler_;
    //! What it was told
    std::vector<ACTRecordedMessage> messages_;
};

/*! \brief A molecule as Python has it
 * With the force field it was made for, which lives as long as the
 * molecule does, and what is needed to compute its energy.
 */
struct Molecule
{
    //! The force field object of Python, which keeps the force field alive
    py::object                     forcefield;
    //! The force field
    ForceField                    *pd = nullptr;
    //! The molecule with its topology and charges
    ACTMol                         mol;
    //! Computes energies and forces
    std::unique_ptr<ForceComputer> forceComp;
    //! The coordinates of all particles in nm
    std::vector<gmx::RVec>         coords;
};

/*! \brief Make a molecule with its topology and charges, as the programs of ACT do
 * \param[in] forcefield The force field object of Python
 * \param[in] fileName   A structure file, empty to make the structure from the SMILES string
 * \param[in] smiles     A SMILES string
 * \return The molecule
 * \throws gmx::InvalidInputError if no molecule can be made of it with this force field
 */
std::unique_ptr<Molecule> makeMolecule(const py::object &forcefield, const std::string &fileName,
                                       const std::string &smiles)
{
    if (!py::isinstance<ForceField>(forcefield))
    {
        throw py::type_error("The force field of a molecule must be an act.ForceField");
    }
    auto made        = std::make_unique<Molecule>();
    made->forcefield = forcefield;
    made->pd         = forcefield.cast<ForceField *>();
    // As alexandria minimize sets it up
    const SimulationConfigHandler sch;
    double                        shellToler = 1e-6;
    if (shellToler >= sch.forceTolerance())
    {
        shellToler = sch.forceTolerance()/10;
    }
    made->forceComp = std::make_unique<ForceComputer>();
    made->forceComp->init(shellToler, sch.shellMaxIter(), made->pd->epsilonr());

    QuietMessages  messages;
    CompoundReader reader;
    const auto     what = fileName.empty() ? "No molecule could be made from the SMILES string '" + smiles + "'" :
        "No molecule could be read from " + fileName;
    if (fileName.empty())
    {
        reader.setSmiles(smiles);
    }
    else
    {
        if (!std::filesystem::exists(fileName))
        {
            GMX_THROW(gmx::FileIOError("Could not find the structure file " + fileName));
        }
        reader.setFileName(fileName);
    }
    std::vector<ACTMol> mols;
    try
    {
        reader.read(messages.handler(), *made->pd, made->forceComp.get(), &mols);
    }
    catch (const gmx::GromacsException &)
    {
        throw;
    }
    catch (const std::exception &ex)
    {
        // The library that reads structures and SMILES strings has exceptions of its own
        GMX_THROW(gmx::InvalidInputError(what + ": " + ex.what()));
    }
    messages.check(what);
    if (mols.empty())
    {
        GMX_THROW(gmx::InvalidInputError(what));
    }
    made->mol    = std::move(mols[0]);
    made->coords = made->mol.xOriginal();
    return made;
}

/*! \brief Compute the energy and the forces of a molecule
 * \param[in]  molecule The molecule
 * \param[out] forces   The force on each particle in kJ/mol nm
 * \return The energy by term in kJ/mol
 * \throws gmx::InvalidInputError if the energy cannot be computed
 */
std::map<InteractionType, double> computeEnergy(const Molecule &molecule, std::vector<gmx::RVec> *forces)
{
    QuietMessages                     messages;
    std::map<InteractionType, double> energies;
    // The shells find their place for these positions of the atoms, in a copy
    auto coords = molecule.coords;
    forces->assign(coords.size(), { 0, 0, 0 });
    molecule.forceComp->compute(messages.handler(), molecule.pd, molecule.mol.topology(), &coords, forces,
                                &energies);
    messages.check("The energy could not be computed");
    return energies;
}

//! \return The coordinates or forces of all particles as a list of lists of three numbers
py::list listOfVectors(const std::vector<gmx::RVec> &vectors)
{
    py::list all;
    for (const auto &vector : vectors)
    {
        py::list one;
        for (int dim = 0; dim < 3; dim++)
        {
            one.append(static_cast<double>(vector[dim]));
        }
        all.append(one);
    }
    return all;
}

//! \return The energies by the names of their terms, for Python
py::dict energyDict(const std::map<InteractionType, double> &energies)
{
    py::dict energy;
    for (const auto &term : energies)
    {
        energy[py::str(interactionTypeToString(term.first))] = term.second;
    }
    return energy;
}

//! The result of a minimization as Python sees it
struct Minimization
{
    //! Whether the forces are within the tolerance
    bool        converged = false;
    //! OK, or why the minimizer stopped, as alexandria minimize tells it
    std::string status;
    //! The distance between the structures before and after in nm, after superposition
    double      rmsd = 0;
    //! The energy by term in kJ/mol before and after
    py::dict    energyBefore, energyAfter;
    //! The position of each particle in nm at the lowest energy found
    py::list    coordinates;
};

/*! \brief Minimize the energy of a molecule, as alexandria minimize does
 * \param[in] molecule  The molecule
 * \param[in] algorithm LBFGS, Steep or Newton, the default of the program if not given
 * \param[in] maxIter   The largest number of iterations, the default of the program if not given
 * \param[in] tolerance The mean square force at which to stop, the default of the program if not given
 * \param[in] freeze    The particles, as in the list of atoms, that stay where they are
 * \param[in] update    Whether the molecule gets the structure that is found
 * \return How it went and the structure
 * \throws gmx::InvalidInputError for settings that make no sense, gmx::ToleranceError
 * if the minimizer failed
 */
Minimization minimizeMolecule(Molecule *molecule, const std::optional<std::string> &algorithm,
                              std::optional<int> maxIter, std::optional<double> tolerance,
                              const std::vector<int> &freeze, bool update)
{
    SimulationConfigHandler sch;
    if (algorithm)
    {
        bool known = false;
        for (auto minAlg : { eMinimizeAlgorithm::LBFGS, eMinimizeAlgorithm::Steep, eMinimizeAlgorithm::Newton })
        {
            if (eMinimizeAlgorithmToString(minAlg) == *algorithm)
            {
                sch.setMinimizeAlgorithm(minAlg);
                known = true;
            }
        }
        if (!known)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown algorithm '%s' to minimize with: use LBFGS, Steep or Newton", algorithm->c_str())));
        }
    }
    if (maxIter)
    {
        if (*maxIter < 1)
        {
            GMX_THROW(gmx::InvalidInputError("The number of iterations of a minimization is at least one"));
        }
        sch.setMaxIter(*maxIter);
    }
    if (tolerance)
    {
        if (!std::isfinite(*tolerance) || *tolerance <= 0)
        {
            GMX_THROW(gmx::InvalidInputError("The tolerance of a minimization is a number above zero"));
        }
        sch.setForceTolerance(*tolerance);
    }
    // The minimizer counts the atoms without shells and virtual sites
    const auto      &atoms = molecule->mol.atomsConst();
    std::vector<int> frozen;
    for (int index : freeze)
    {
        if (index < 0 || index >= static_cast<int>(atoms.size()))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("There is no particle %d to freeze, the molecule has %zu", index, atoms.size())));
        }
        if (atoms[index].pType() != ActParticle::Atom)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Particle %d is a %s, only atoms can be frozen", index,
                                                               actParticleToString(atoms[index].pType()).c_str())));
        }
        int number = 0;
        for (int i = 0; i < index; i++)
        {
            number += (atoms[i].pType() == ActParticle::Atom) ? 1 : 0;
        }
        frozen.push_back(number);
    }
    // The shells more precisely than the atoms, as the program has it
    double shellToler = 1e-6;
    if (shellToler >= sch.forceTolerance())
    {
        shellToler = sch.forceTolerance()/10;
    }
    molecule->forceComp->init(shellToler, sch.shellMaxIter(), molecule->pd->epsilonr());

    Minimization           result;
    std::vector<gmx::RVec> forces;
    result.energyBefore = energyDict(computeEnergy(*molecule, &forces));

    QuietMessages                     messages;
    std::map<InteractionType, double> energies;
    auto                              xmin = molecule->coords;
    const MolHandler                  molhandler;
    const auto status = molhandler.minimizeCoordinates(messages.handler(), molecule->pd, molecule->mol.topology(),
                                                       molecule->forceComp.get(), sch, &xmin, &energies, frozen);
    messages.check("The minimization failed");
    if (eMinimizeStatus::OK != status && !minimizationDidNotConverge(status))
    {
        GMX_THROW(gmx::ToleranceError("The minimization failed: " +
                                      minimizeStatusMessage(status, energies[InteractionType::EPOT])));
    }
    result.converged   = (eMinimizeStatus::OK == status);
    result.status      = eMinimizeStatusToString(status);
    // On top of the structure it started from
    result.rmsd        = molhandler.coordinateRmsd(atoms, molecule->coords, &xmin);
    result.energyAfter = energyDict(energies);
    result.coordinates = listOfVectors(xmin);
    if (update)
    {
        molecule->coords = xmin;
    }
    return result;
}

//! Tell Python about the result of a minimization
void defineMinimization(py::module_ &module)
{
    py::class_<Minimization>(module, "Minimization", "The result of Molecule.minimize.")
        .def_readonly("converged", &Minimization::converged,
                      "Whether the forces are within the tolerance.")
        .def_readonly("status", &Minimization::status,
                      "OK, or why the minimizer stopped, such as Too many steps. The structure is\n"
                      "then the one with the lowest energy that was found.")
        .def_readonly("rmsd", &Minimization::rmsd,
                      "The root mean square distance in nm between the atoms before and after,\n"
                      "with the structures on top of each other.")
        .def_readonly("energy_before", &Minimization::energyBefore,
                      "The energy by term in kJ/mol of the structure it started from.")
        .def_readonly("energy_after", &Minimization::energyAfter,
                      "The energy by term in kJ/mol of the structure that was found.")
        .def_readonly("coordinates", &Minimization::coordinates,
                      "The position of each particle in nm in the structure that was found.")
        .def("__repr__",
             [](const Minimization &result)
             {
                 return gmx::formatString("<act.Minimization %s (%s), RMSD %g nm>",
                                          result.converged ? "converged" : "not converged",
                                          result.status.c_str(), result.rmsd);
             });
}

//! The vibrations of a molecule as Python sees them
struct NormalModes
{
    //! The frequencies in cm^-1, from low to high, and the infrared intensities
    std::vector<double> frequencies, intensities;
    //! The unit of the intensities
    std::string         intensityUnit;
    //! For each mode its frequency, intensity, reduced mass and the displacement of each particle
    py::list            modes;
    //! Zero point energy, entropy, heat capacity and so on, as alexandria nma reports them
    py::object          thermochemistry;
    //! The minimization that was done first, None if there was none
    py::object          minimization;
};

/*! \brief The normal modes of a molecule, as alexandria nma computes them
 * \param[in] molecule  The molecule
 * \param[in] minimize  Whether to minimize the energy first, which the molecule keeps
 * \param[in] algorithm For the minimization, see minimizeMolecule
 * \param[in] maxIter   For the minimization
 * \param[in] tolerance For the minimization
 * \return The vibrations
 * \throws gmx::ToleranceError if the minimization does not converge,
 * gmx::InvalidInputError if the analysis fails
 */
NormalModes normalModes(Molecule *molecule, bool minimize, const std::optional<std::string> &algorithm,
                        std::optional<int> maxIter, std::optional<double> tolerance)
{
    NormalModes result;
    result.minimization = py::none();
    if (minimize)
    {
        const auto minimum = minimizeMolecule(molecule, algorithm, maxIter, tolerance, {}, true);
        if (!minimum.converged)
        {
            GMX_THROW(gmx::ToleranceError("The minimization before the normal modes did not converge (" + minimum.status +
                                          "), and vibrations of a structure that is no minimum make no sense. "
                                          "Allow more iterations, or minimize the molecule first"));
        }
        result.minimization = py::cast(minimum);
    }
    QuietMessages     messages;
    AtomizationEnergy atomenergy;
    atomenergy.read(messages.handler());
    messages.check("The atomization energies of ACT could not be read");
    JsonTree         jtree("nma");
    auto             coords = molecule->coords;
    const MolHandler molhandler;
    // Without a file for the spectrum, and with the motion of the particles
    doFrequencyAnalysis(messages.handler(), molecule->pd, &molecule->mol, molhandler, molecule->forceComp.get(),
                        &coords, atomenergy, nullptr, &jtree, nullptr, 24, nullptr, false, true);
    messages.check("The normal modes could not be computed");
    const py::dict told = toPython(jtree);
    if (!told.contains("modes"))
    {
        GMX_THROW(gmx::InternalError("The normal mode analysis gave no modes"));
    }
    result.modes = py::list();
    for (const auto &item : told["modes"])
    {
        const auto mode = py::reinterpret_borrow<py::dict>(item);
        py::dict   one;
        one["frequency"]     = mode["frequency"]["value"];
        one["intensity"]     = mode["intensity"]["value"];
        one["reduced_mass"]  = mode["reduced_mass"]["value"];
        one["displacements"] = mode["displacements"];
        result.modes.append(one);
        result.frequencies.push_back(mode["frequency"]["value"].cast<double>());
        result.intensities.push_back(mode["intensity"]["value"].cast<double>());
        result.intensityUnit = mode["intensity"]["unit"].cast<std::string>();
    }
    result.thermochemistry = told["Thermochemistry"]["Alexandria"];
    return result;
}

/*! \brief The infrared spectrum of the vibrations, as alexandria nma -ir writes it
 * A Lorentzian at each frequency, for whole numbers of cm^-1 from five line
 * widths below the lowest to five above the highest frequency, and
 * normalized such that the intensities add up to one.
 * \param[in] modes     The vibrations
 * \param[in] linewidth The width of the lines in cm^-1
 * \return The frequencies in cm^-1 and the intensities
 * \throws gmx::InvalidInputError if the line width is not above zero
 */
py::tuple spectrumOf(const NormalModes &modes, double linewidth)
{
    if (!std::isfinite(linewidth) || linewidth <= 0)
    {
        GMX_THROW(gmx::InvalidInputError("The line width of a spectrum is a number above zero"));
    }
    std::vector<int>    frequencies;
    std::vector<double> spectrum;
    if (!modes.frequencies.empty())
    {
        const int minf = *std::min_element(modes.frequencies.begin(), modes.frequencies.end()) - 5*linewidth;
        const int maxf = *std::max_element(modes.frequencies.begin(), modes.frequencies.end()) + 5*linewidth;
        spectrum.assign(maxf - minf + 1, 0.0);
        double total = 0;
        for (int f = minf; f <= maxf; f++)
        {
            frequencies.push_back(f);
            for (size_t peak = 0; peak < modes.frequencies.size(); peak++)
            {
                spectrum[f - minf] += modes.intensities[peak]*linewidth/
                    (2*M_PI*(gmx::square(f - modes.frequencies[peak]) + gmx::square(linewidth/2)));
            }
            total += spectrum[f - minf];
        }
        if (total > 0)
        {
            for (auto &value : spectrum)
            {
                value /= total;
            }
        }
    }
    return py::make_tuple(frequencies, spectrum);
}

//! Tell Python about the vibrations of a molecule
void defineNormalModes(py::module_ &module)
{
    py::class_<NormalModes>(module, "NormalModes", "The result of Molecule.normal_modes.")
        .def_readonly("frequencies", &NormalModes::frequencies,
                      "The frequencies of the vibrations in cm^-1, from low to high.")
        .def_readonly("intensities", &NormalModes::intensities,
                      "The infrared intensity of each vibration, in intensity_unit.")
        .def_readonly("intensity_unit", &NormalModes::intensityUnit, "The unit of the intensities.")
        .def_readonly("modes", &NormalModes::modes,
                      "For each vibration its frequency, intensity, reduced_mass in g/mol and\n"
                      "displacements: how each particle moves, as x, y and z, normalized as\n"
                      "quantum chemistry programs print them.")
        .def_readonly("thermochemistry", &NormalModes::thermochemistry,
                      "The zero point energy, entropy, heat capacity and so on, each with value\n"
                      "and unit, as alexandria nma reports them.")
        .def_readonly("minimization", &NormalModes::minimization,
                      "The Minimization that was done first, None if there was none.")
        .def("spectrum", &spectrumOf, py::arg("linewidth") = 24.0,
             "The infrared spectrum: a line of the given width in cm^-1 at each frequency.\n\n"
             "Returns the frequencies, whole numbers of cm^-1, and the intensities, which\n"
             "add up to one, as alexandria nma -ir writes them.")
        .def("__repr__",
             [](const NormalModes &modes)
             {
                 return gmx::formatString("<act.NormalModes %zu modes, %g to %g cm^-1>", modes.frequencies.size(),
                                          modes.frequencies.empty() ? 0.0 : modes.frequencies.front(),
                                          modes.frequencies.empty() ? 0.0 : modes.frequencies.back());
             });
}

//! A simulation as Python sees it
struct Trajectory
{
    //! The number of steps that were asked for and the time step in ps
    int                                        steps = 0;
    double                                     timestep = 0;
    //! The steps and times in ps of the energies
    std::vector<int>                           energySteps;
    std::vector<double>                        energyTimes;
    //! The temperature in K at these times
    std::vector<double>                        temperature;
    //! The energies in kJ/mol at these times, by term, with EKIN and ETOT
    std::map<std::string, std::vector<double>> energies;
    //! The steps and times in ps of the structures
    std::vector<int>                           frameSteps;
    std::vector<double>                        times;
    //! The positions of all particles in nm at these times
    std::vector<std::vector<gmx::RVec>>        frames;
};

/*! \brief Run a molecular dynamics simulation of a molecule, as alexandria simulate does
 * \param[in] molecule       The molecule, of which the structure stays as it is
 * \param[in] steps          The number of steps
 * \param[in] timestep       The time step in ps
 * \param[in] temperature    The temperature in K of the velocities at the start, zero to start at rest
 * \param[in] seed           The seed for the velocities, zero to make one
 * \param[in] energyInterval The number of steps between energies, zero for none
 * \param[in] frameInterval  The number of steps between structures, zero for none
 * \param[in] onEnergies     A function of Python that gets the energies when they are due
 *                           and ends the simulation by returning False, or None
 * \return The energies and structures along the way
 * \throws gmx::InvalidInputError for settings that make no sense
 */
Trajectory simulateMolecule(Molecule &molecule, int steps, double timestep, double temperature, int seed,
                            int energyInterval, int frameInterval, const py::object &onEnergies)
{
    if (steps < 1)
    {
        GMX_THROW(gmx::InvalidInputError("A simulation has at least one step"));
    }
    if (!std::isfinite(timestep) || timestep <= 0)
    {
        GMX_THROW(gmx::InvalidInputError("The time step of a simulation is a number above zero, in ps"));
    }
    if (!std::isfinite(temperature) || temperature < 0)
    {
        GMX_THROW(gmx::InvalidInputError("The temperature of a simulation cannot be negative"));
    }
    if (energyInterval < 0 || frameInterval < 0 || seed < 0)
    {
        GMX_THROW(gmx::InvalidInputError("The seed and the numbers of steps between energies and between structures cannot be negative"));
    }
    if (!onEnergies.is_none() && !PyCallable_Check(onEnergies.ptr()))
    {
        throw py::type_error("on_energies must be a function, or None");
    }
    SimulationConfigHandler sch;
    sch.setSimulation(steps, timestep, temperature, seed, energyInterval, frameInterval);

    Trajectory result;
    result.steps    = steps;
    result.timestep = timestep;
    // An error of the function of Python, or Ctrl-C, ends the simulation
    // and is raised when the simulation has returned
    bool               failed = false;
    SimulationObserver observer;
    observer.energies = [&](int step, double time, double temp, double ekin,
                            const std::map<InteractionType, double> &energies)
    {
        result.energySteps.push_back(step);
        result.energyTimes.push_back(time);
        result.temperature.push_back(temp);
        double epot = 0;
        for (const auto &term : energies)
        {
            result.energies[interactionTypeToString(term.first)].push_back(term.second);
            if (InteractionType::EPOT == term.first)
            {
                epot = term.second;
            }
        }
        result.energies["EKIN"].push_back(ekin);
        result.energies["ETOT"].push_back(ekin + epot);
        const py::gil_scoped_acquire python;
        if (PyErr_CheckSignals() != 0)
        {
            failed = true;
            return false;
        }
        if (onEnergies.is_none())
        {
            return true;
        }
        try
        {
            py::dict item;
            item["step"]        = step;
            item["time"]        = time;
            item["temperature"] = temp;
            py::dict values;
            for (const auto &term : result.energies)
            {
                values[py::str(term.first)] = term.second.back();
            }
            item["energies"] = values;
            const py::object answer = onEnergies(item);
            // Only False ends the simulation: a function that returns nothing lets it run
            return !(py::isinstance<py::bool_>(answer) && !answer.cast<bool>());
        }
        catch (py::error_already_set &error)
        {
            error.restore();
            failed = true;
            return false;
        }
    };
    observer.coordinates = [&](int step, double time, const std::vector<gmx::RVec> &coords)
    {
        result.frameSteps.push_back(step);
        result.times.push_back(time);
        result.frames.push_back(coords);
        return true;
    };
    QuietMessages    messages;
    const MolHandler molhandler;
    {
        // Other threads of Python run while this one computes. The
        // simulation starts from the coordinates that are given and keeps
        // its own, so the structure of the molecule is not changed.
        const py::gil_scoped_release release;
        molhandler.simulate(messages.handler(), molecule.pd, &molecule.mol, molecule.coords,
                            molecule.forceComp.get(), sch, nullptr, nullptr, nullptr, nullptr, &observer);
    }
    if (failed)
    {
        throw py::error_already_set();
    }
    messages.check("The simulation failed");
    return result;
}

//! Tell Python about the result of a simulation
void defineTrajectory(py::module_ &module)
{
    py::class_<Trajectory>(module, "Trajectory", "The result of Molecule.simulate.")
        .def_readonly("steps", &Trajectory::steps, "The number of steps that were asked for.")
        .def_readonly("timestep", &Trajectory::timestep, "The time step in ps.")
        .def_readonly("energy_steps", &Trajectory::energySteps, "The steps of which there are energies.")
        .def_readonly("energy_times", &Trajectory::energyTimes, "The times in ps of the energies.")
        .def_readonly("temperature", &Trajectory::temperature, "The temperature in K at these times.")
        .def_readonly("energies", &Trajectory::energies,
                      "The energies in kJ/mol at these times: for each term of the force field, for\n"
                      "the potential (EPOT), kinetic (EKIN) and total energy (ETOT) a list.")
        .def_readonly("frame_steps", &Trajectory::frameSteps, "The steps of which there is a structure.")
        .def_readonly("times", &Trajectory::times, "The times in ps of the structures.")
        .def_property_readonly("frames",
                               [](const Trajectory &trajectory)
                               {
                                   py::list frames;
                                   for (const auto &frame : trajectory.frames)
                                   {
                                       frames.append(listOfVectors(frame));
                                   }
                                   return frames;
                               },
                               "The structures: for each the position of every particle in nm.")
        .def("__repr__",
             [](const Trajectory &trajectory)
             {
                 return gmx::formatString("<act.Trajectory %d steps of %g ps, %zu energies, %zu structures>",
                                          trajectory.steps, trajectory.timestep, trajectory.energySteps.size(),
                                          trajectory.frames.size());
             });
}

/*! \brief Tell Python about the class Molecule
 * A lambda that returns one item of what toPython gives has to say that it
 * returns an object: an item of a temporary is otherwise looked up after
 * the temporary is gone.
 */
void defineMolecule(py::module_ &module)
{
    py::class_<Molecule>(module, "Molecule",
                         "A molecule for a force field: its particles with their charges, its bonds,\n"
                         "its structure and its energy.\n\n"
                         "Make one with Molecule.from_smiles(smiles, forcefield) or\n"
                         "Molecule.from_file(file, forcefield).")
        .def_static("from_smiles",
                    [](const std::string &smiles, const py::object &forcefield)
                    {
                        return makeMolecule(forcefield, "", smiles);
                    },
                    py::arg("smiles"), py::arg("forcefield"),
                    "Make a molecule from a SMILES string: a structure, the topology for the\n"
                    "force field and the charges, as the programs of ACT make them.\n\n"
                    "Raises ActInputError if the string is no molecule or the force field has\n"
                    "no types for its atoms.")
        .def_static("from_file",
                    [](const std::filesystem::path &fileName, const py::object &forcefield)
                    {
                        return makeMolecule(forcefield, fileName.string(), "");
                    },
                    py::arg("file"), py::arg("forcefield"),
                    "Read a molecule from a structure file, such as an sdf or pdb file, and make\n"
                    "the topology for the force field and the charges.\n\n"
                    "Raises ActFileError if the file is not there, and ActInputError if no\n"
                    "molecule can be made of it with this force field.")
        .def_readonly("forcefield", &Molecule::forcefield, "The force field the molecule was made for.")
        .def_property_readonly("formula",
                               [](const Molecule &molecule) -> py::object
                               {
                                   return toPython(moleculeToJson(&molecule.mol))["formula"];
                               },
                               "The chemical formula.")
        .def_property_readonly("name",
                               [](const Molecule &molecule) -> py::object
                               {
                                   return toPython(moleculeToJson(&molecule.mol))["name"];
                               },
                               "The name of the molecule, empty if it has none.")
        .def_property_readonly("charge",
                               [](const Molecule &molecule) -> py::object
                               {
                                   return toPython(moleculeToJson(&molecule.mol))["charge"];
                               },
                               "The total charge in e.")
        .def_property_readonly("mass",
                               [](const Molecule &molecule) -> py::object
                               {
                                   return toPython(moleculeToJson(&molecule.mol))["mass"]["value"];
                               },
                               "The mass in g/mol.")
        .def_property_readonly("atoms",
                               [](const Molecule &molecule) -> py::object
                               {
                                   return toPython(moleculeToJson(&molecule.mol))["particles"];
                               },
                               "The particles: atoms, shells and virtual sites. For each its name,\n"
                               "element, type in the force field, kind (particletype), mass and charge,\n"
                               "in the order of the coordinates.")
        .def_property_readonly("bonds",
                               [](const Molecule &molecule) -> py::object
                               {
                                   return toPython(moleculeToJson(&molecule.mol))["bonds"];
                               },
                               "The bonds: for each the two particles, counted from zero, and the bond order.")
        .def_property_readonly("charges",
                               [](const Molecule &molecule)
                               {
                                   py::list charges;
                                   for (const auto &atom : molecule.mol.atomsConst())
                                   {
                                       charges.append(static_cast<double>(atom.charge()));
                                   }
                                   return charges;
                               },
                               "The charge of each particle in e.")
        .def_property_readonly("coordinates",
                               [](const Molecule &molecule) { return listOfVectors(molecule.coords); },
                               "The position of each particle in nm, as a list of x, y and z.")
        .def("energy",
             [](const Molecule &molecule)
             {
                 std::vector<gmx::RVec> forces;
                 return energyDict(computeEnergy(molecule, &forces));
             },
             "The energy in kJ/mol for the structure of the molecule: a dict with the terms\n"
             "of the force field and their sum, EPOT.")
        .def("forces",
             [](const Molecule &molecule)
             {
                 std::vector<gmx::RVec> forces;
                 (void) computeEnergy(molecule, &forces);
                 return listOfVectors(forces);
             },
             "The force on each particle in kJ/mol nm, as a list of x, y and z.")
        .def("minimize", &minimizeMolecule, py::kw_only(), py::arg("algorithm") = py::none(),
             py::arg("max_iter") = py::none(), py::arg("tolerance") = py::none(),
             py::arg("freeze") = std::vector<int>(), py::arg("update") = true,
             "Minimize the energy with respect to the positions of the atoms, as alexandria\n"
             "minimize does, and give the molecule the structure that is found.\n\n"
             "algorithm is LBFGS, Steep or Newton, max_iter the largest number of iterations\n"
             "and tolerance the mean square force at which to stop; what is not given is as\n"
             "in the program. The particles in freeze, numbered as in atoms, stay where they\n"
             "are. With update=False the molecule keeps its structure. Returns a Minimization\n"
             "with whether it converged, the energies before and after, the RMSD and the\n"
             "structure. A minimization that does not converge is told there, not raised.")
        .def("normal_modes", &normalModes, py::kw_only(), py::arg("minimize") = true,
             py::arg("algorithm") = py::none(), py::arg("max_iter") = py::none(),
             py::arg("tolerance") = py::none(),
             "The vibrations of the molecule, as alexandria nma computes them: frequencies,\n"
             "infrared intensities, the motion of the particles and the thermochemistry.\n\n"
             "The energy is minimized first, with algorithm, max_iter and tolerance as in\n"
             "minimize, and the molecule keeps that structure; give minimize=False for a\n"
             "molecule that is at its minimum already. Raises ActError if the minimization\n"
             "does not converge. Returns a NormalModes.")
        .def("simulate", &simulateMolecule, py::kw_only(), py::arg("steps"), py::arg("timestep") = 0.0002,
             py::arg("temperature") = 0.0, py::arg("seed") = 0, py::arg("energy_interval") = 1,
             py::arg("frame_interval") = 1, py::arg("on_energies") = py::none(),
             "Run a molecular dynamics simulation, as alexandria simulate does, from the\n"
             "structure of the molecule, which stays as it is.\n\n"
             "steps is the number of steps of timestep ps. With a temperature in K the atoms\n"
             "get velocities for it, drawn with the seed (0 for another one each time);\n"
             "without, they start at rest. The energies are kept every energy_interval steps\n"
             "and the structure every frame_interval steps, 0 for never. on_energies is a\n"
             "function that gets a dict with step, time, temperature and energies each time\n"
             "the energies are due; if it returns False the simulation ends there. Returns a\n"
             "Trajectory. Ctrl-C ends the simulation at the next energies.")
        .def("__repr__",
             [](const Molecule &molecule)
             {
                 const auto described = toPython(moleculeToJson(&molecule.mol));
                 return "<act.Molecule " + described["formula"].cast<std::string>() + ", " +
                        std::to_string(molecule.mol.atomsConst().size()) + " particles>";
             });
}

} // namespace

} // namespace alexandria

PYBIND11_MODULE(_core, module)
{
    module.doc() = "The classes of ACT in Python. Use them through the package act.";
    // A library must leave its caller alive: errors are exceptions from here on
    gmx_set_fatal_error_throws(true);
    // The data files of ACT are where this module is installed
    gmx::setProgramContext(&alexandria::moduleContext());
    module.def("data_directory", []() { return alexandria::moduleContext().dataDirectory(); },
               "The directory of the installation in which ACT looks for its data files, such as\n"
               "the properties of the elements, after the directory in the environment variable\n"
               "ACTDATA. Empty if the installation has none.");
    py::register_exception_translator(&alexandria::translateException);
    alexandria::defineParameters(module);
    alexandria::defineForceField(module);
    alexandria::defineMinimization(module);
    alexandria::defineNormalModes(module);
    alexandria::defineTrajectory(module);
    alexandria::defineMolecule(module);
}
