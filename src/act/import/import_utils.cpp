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
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
    
#include "import_utils.h"

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "act/basics/msg_handler.h"
#include "act/utility/xml_util.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

//! \brief Enum for distinguishing content types in molprop files.
enum class AtomBondtypeXml {
    ABTYPE,
    NAME,
    SMARTS,
    CHARGE,
    MULTIPLICITY,
    INDEX,
    ATOMTYPES,
    ATOMTYPE,
    BONDTYPE,
    BONDTYPES,
    AI,
    AJ,
    ORDER
};

//! \brief Map string to AtomBondtypeXml
static std::map<const std::string, AtomBondtypeXml> xmlxxx =
{
    { "abtype", AtomBondtypeXml::ABTYPE },
    { "name", AtomBondtypeXml::NAME },
    { "smarts", AtomBondtypeXml::SMARTS },
    { "charge", AtomBondtypeXml::CHARGE },
    { "multiplicity", AtomBondtypeXml::MULTIPLICITY },
    { "index", AtomBondtypeXml::INDEX },
    { "atomtypes", AtomBondtypeXml::ATOMTYPES },
    { "atomtype", AtomBondtypeXml::ATOMTYPE },
    { "bondtype", AtomBondtypeXml::BONDTYPE },
    { "bondtypes", AtomBondtypeXml::BONDTYPES },
    { "ai", AtomBondtypeXml::AI },
    { "aj", AtomBondtypeXml::AJ },
    { "order", AtomBondtypeXml::ORDER }
};

//! Map MolPropXml to string
static std::map<AtomBondtypeXml, const std::string> rmap = {};

//! \brief Utility to fill utility structures
static void fillMaps()
{
    rmap.clear();
    for (auto &xo : xmlxxx)
    {
        rmap.insert({ xo.second, xo.first });
    }
}

static void mp_process_tree(gmx_unused MsgHandler                             *msghandler,
                            gmx_unused xmlNodePtr                              parent,
                            gmx_unused std::vector<AtomBondtypeEntry>         *abdb,
                            gmx_unused std::map<AtomBondtypeXml, std::string> *xbuf)
{
}

void readAtomBondtypeDB(MsgHandler                     *msghandler,
                        const std::string              &dbname,
                        std::vector<AtomBondtypeEntry> *abdb)
{
    xmlDocPtr   doc;
    std::string mpfile;

    fillMaps();
    if (msghandler)
    {
        msghandler->msg(ACTStatus::Verbose, gmx::formatString("Opening %s\n", dbname.c_str()));
    }

    if ((doc = xmlParseFile(dbname.c_str())) == nullptr)
    {
        if (msghandler)
        {
            msghandler->msg(ACTStatus::Error,
                            gmx::formatString("Failed reading XML file %s. Run a syntax checker such as nsgmls.",
                                              dbname.c_str()));
        }
    }

    std::map<AtomBondtypeXml, std::string>  xbuf;
    mp_process_tree(msghandler, doc->children, abdb, &xbuf);
    xmlFreeDoc(doc);
}

std::vector<AtomBondtypeEntry> getAtomBondtypeDB()
{
    // Note that the order is important, the vector will be checked from top to bottom
    std::vector<AtomBondtypeEntry> abe = {
        { "sulfate",
          "[#16](=[#8])(=[#8])(-[#8-])-[#8-]",
          -2, 1,
          { "s3", "o2", "o2", "o2", "o2" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 }, { 0, 3, 1.5 }, { 0, 4, 1.5 } }
        },
        { "water",
          "[#8](-[#1])-[#1]",
          0, 1,
          { "ow", "hw", "hw" },
          { { 0, 1, 1 }, { 0, 2, 1 } }
        },
        { "phosphate", // hypervalent form
          "[#15D4]([#8D1])([#8D1])([#8-,#8D1])([#8-,#8D1])",
          -3, 1,
          { "p3", "o3", "o3", "o3", "o3" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 }, { 0, 3, 1.5 }, { 0, 4, 1.5 } }
        },
        { "phosphate2", // ion form, PO4(3-)
          "[#8-]-[#15](-[#8-])(-[#8-])=[#8]",
          -3, 1,
          { "o3", "p3", "o3", "o3", "o3" },
          { { 0, 1, 1.5 }, { 1, 2, 1.5 }, { 1, 3, 1.5 }, { 1, 4, 1.5 } }
        },
        { "phosphate2", // ion form, XPO3(2-)
          "[#8-]-[#15](-[#8-])(-[#8-])",
          -2, 1,
          { "o3", "p3", "o3", "o3" },
          { { 0, 1, 1.5 }, { 1, 2, 1.5 }, { 1, 3, 1.5 } }
        },
        { "phosphate3", // ion form, X2PO2(1-)
          "[#8-]-[#15](-[#8-])",
          -1, 1,
          { "o3", "p3", "o3" },
          { { 0, 1, 1 }, { 1, 2, 1 } }
        },
        { "nitro1",
          "[#7+](-[#8-])=[#8]",
          0, 1,
          { "n2", "o2", "o2" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 } }
        },
        { "nitro2",
          "[#7+](-[#8])=[#8]",
          0, 1,
          { "n2", "o2", "o2" },
          { { 0, 1, 1.5 }, { 0, 2, 1.5 } }
        }
    };
    return abe;
}

static void add_xml_abentry(xmlNodePtr               parent,
                            const AtomBondtypeEntry &ab)
{
    xmlNodePtr ptr = add_xml_child(parent, rmap[AtomBondtypeXml::ABTYPE]);
    add_xml_string(ptr, rmap[AtomBondtypeXml::NAME], ab.name);
    add_xml_string(ptr, rmap[AtomBondtypeXml::SMARTS], ab.smarts);
    add_xml_int(ptr, rmap[AtomBondtypeXml::CHARGE], ab.charge);
    add_xml_int(ptr, rmap[AtomBondtypeXml::MULTIPLICITY], ab.multiplicity);
    auto atypes = add_xml_child(ptr, rmap[AtomBondtypeXml::ATOMTYPES]);
    int index = 0;
    for(const auto &atp : ab.atomtypes)
    {
        auto myatp = add_xml_child(atypes, rmap[AtomBondtypeXml::ATOMTYPE]);
        add_xml_int(myatp, rmap[AtomBondtypeXml::INDEX], index);
        add_xml_string(myatp, rmap[AtomBondtypeXml::NAME], atp.c_str());
        index += 1;
    }
    auto btypes = add_xml_child(ptr, rmap[AtomBondtypeXml::BONDTYPES]);
    for(const auto &btp : ab.bonds)
    {
        auto mybtp = add_xml_child(btypes, rmap[AtomBondtypeXml::BONDTYPE]);
        add_xml_int(mybtp, rmap[AtomBondtypeXml::AI], btp.ai);
        add_xml_int(mybtp, rmap[AtomBondtypeXml::AJ], btp.aj);
        add_xml_double(mybtp, rmap[AtomBondtypeXml::ORDER], btp.order);
    }
}

void writeAtomBondtypeDB(MsgHandler                                      *msghandler,
                         gmx_unused const std::string                    &filenm,
                         gmx_unused const std::vector<AtomBondtypeEntry> &db)
{
    xmlDocPtr   doc;
    xmlDtdPtr   dtd;
    xmlNodePtr  myroot;
    xmlChar    *libdtdname, *dtdname, *gmx;

    gmx        = (xmlChar *) "AtomBondtype";
    dtdname    = (xmlChar *) "atombondDB.dtd";
    libdtdname = dtdname;
    fillMaps();
    if ((doc = xmlNewDoc((xmlChar *)"1.0")) == nullptr)
    {
        msghandler->msg(ACTStatus::Error, gmx::formatString("Creating XML document %s", filenm.c_str()));
    }

    if ((dtd = xmlCreateIntSubset(doc, dtdname, libdtdname, dtdname)) == nullptr)
    {
        msghandler->msg(ACTStatus::Error, gmx::formatString("Creating XML DTD %s", filenm.c_str()));
    }

    if ((myroot = xmlNewDocNode(doc, nullptr, gmx, nullptr)) == nullptr)
    {
        msghandler->msg(ACTStatus::Error, gmx::formatString("Creating root element for %s", filenm.c_str()));
    }
    dtd->next    = myroot;
    myroot->prev = (xmlNodePtr) dtd;
    // Recreate reverse map
    fillMaps();

    // Add AtomBondtype definitions
    for (auto &ab : db)
    {
        add_xml_abentry(myroot, ab);
    }
    // No compression
    xmlSetDocCompressMode(doc, 0);
    xmlIndentTreeOutput = 1;
    if (xmlSaveFormatFileEnc(filenm.c_str(), doc, "ISO-8859-1", 1) == 0)
    {
        msghandler->msg(ACTStatus::Error, gmx::formatString("Faild to save file %s", filenm.c_str()));
    }
    xmlFreeDoc(doc);
}

} // namespace alexandria
