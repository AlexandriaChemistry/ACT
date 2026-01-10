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
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

//! \brief Enum for distinguishing content types in molprop files.
enum class AtomBondtypeXml {
    ABTYPES,
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
    { "atombondtypes", AtomBondtypeXml::ABTYPES },
    { "atombondtype", AtomBondtypeXml::ABTYPE },
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

//! Map AtomBondtypeXml to string
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

//! \brief Clean specific elements from XML buffer
static void clean_xbuf(std::map<AtomBondtypeXml, std::string> *xbuf,
                       const std::vector<AtomBondtypeXml>     &clean)
{
    for (auto &mpx : clean)
    {
        auto ff = xbuf->find(mpx);
        if (xbuf->end() != ff)
        {
            xbuf->erase(ff);
        }
    }
}

//! \brief Extract double from xml buffer map or crash
static double xbuf_atof(const std::map<AtomBondtypeXml, std::string> *xbuf, AtomBondtypeXml index, bool crash)
{
    auto ptr = xbuf->find(index);
    if (ptr != xbuf->end())
    {
        return my_atof(ptr->second.c_str(), rmap[index].c_str());
    }
    if (crash)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Input for %s is missing", 
                                                           rmap[index].c_str()).c_str()));
    }
    return 0.0;
}

//! \brief Extract int from xml buffer map or crash
static int xbuf_atoi(const std::map<AtomBondtypeXml, std::string> *xbuf, AtomBondtypeXml index, bool crash)
{
    auto ptr = xbuf->find(index);
    if (ptr != xbuf->end())
    {
        return my_atoi(ptr->second.c_str(), rmap[index].c_str());
    }
    if (crash)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Input for %s is missing", 
                                                           rmap[index].c_str()).c_str()));
    }
    return 0;
}

//! \return true if all AtomBondtypeXml in the index have been read
static bool xmlFound(const std::map<AtomBondtypeXml, std::string> *xbuf, 
                     const std::vector<AtomBondtypeXml>           &index)
{
    bool found = true;
    for(auto &ind : index)
    {
        auto ptr = xbuf->find(ind);
        found = found && (xbuf->end() != ptr) && !ptr->second.empty();
    }
    return found;
}

//! \brief Get XML attributes for an entry
static void get_attributes(MsgHandler                             *msg_handler,
                           xmlAttrPtr                              attr,
                           std::map<AtomBondtypeXml, std::string> *xbuf)
{
    while (attr != nullptr)
    {
        std::string attrname(reinterpret_cast<const char *>(attr->name));
        std::string attrval(reinterpret_cast<const char *>(attr->children->content));

        auto  iter = xmlxxx.find(attrname);
        if (iter != xmlxxx.end() && !attrval.empty())
        {
            auto ptr = xbuf->find(iter->second);
            if (xbuf->end() == ptr)
            {
                xbuf->insert({iter->second, attrval});
            }
            else
            {
                ptr->second.assign(attrval);
            }
        }
        if (msg_handler)
        {
            msg_handler->writeDebug(gmx::formatString("Property: '%s' Value: '%s'\n",
                                                      attrname.c_str(), attrval.c_str()));
        }
        attr = attr->next;
    }
}

static void mp_process_tree(MsgHandler                             *msghandler,
                            xmlNodePtr                              tree,
                            std::vector<AtomBondtypeEntry>         *abdb,
                            std::map<AtomBondtypeXml, std::string> *xbuf)
{
    auto xmltype = xmltypes();
    while (tree != nullptr)
    {
        if (msghandler)
        {
            if ((tree->type > 0) && ((unsigned)tree->type < xmltype.size()))
            {
                msghandler->msg(ACTStatus::Debug,
                                 gmx::formatString("Node type %s encountered with name %s\n",
                                                   xmltype[tree->type], (char *)tree->name));
            }
            else
            {
                msghandler->msg(ACTStatus::Warning,
                                 gmx::formatString("Unknown node type %d encountered\n", tree->type));
            }
        }
        if (tree->type == XML_ELEMENT_NODE)
        {
            auto iter = xmlxxx.find((const char *)tree->name);
            if (iter != xmlxxx.end())
            {
                AtomBondtypeXml elem = iter->second;
                if (msghandler)
                {
                    msghandler->writeDebug(gmx::formatString("Element node name %s\n", (char *)tree->name));
                }
                if (tree->children && tree->children->content)
                {
                    xbuf->insert({ elem, (const char *)tree->children->content });
                }
                get_attributes(msghandler, tree->properties, xbuf);
                
                switch (elem)
                {
                case AtomBondtypeXml::ABTYPES:
                    mp_process_tree(msghandler, tree->children, abdb, xbuf);
                    clean_xbuf(xbuf, { elem });
                    break;
                case AtomBondtypeXml::ABTYPE:
                    {
                        std::vector<AtomBondtypeXml> clean = {
                            AtomBondtypeXml::NAME, AtomBondtypeXml::SMARTS,
                            AtomBondtypeXml::CHARGE, AtomBondtypeXml::MULTIPLICITY
                        };
                        if (xmlFound(xbuf, clean))
                        {
                            AtomBondtypeEntry abe;
                            abe.name.assign(xbuf->find(AtomBondtypeXml::NAME)->second);
                            abe.smarts.assign(xbuf->find(AtomBondtypeXml::SMARTS)->second);
                            abe.charge = xbuf_atoi(xbuf, AtomBondtypeXml::CHARGE, true);
                            abe.multiplicity = xbuf_atoi(xbuf, AtomBondtypeXml::MULTIPLICITY, true);
                            abdb->push_back(std::move(abe));
                            clean_xbuf(xbuf, clean);
                        }
                        else
                        {
                            GMX_THROW(gmx::InvalidInputError("Incomplete AtomBondtype definition"));
                        }
                    }
                    mp_process_tree(msghandler, tree->children, abdb, xbuf);
                    break;
                case AtomBondtypeXml::ATOMTYPES:
                case AtomBondtypeXml::BONDTYPES:
                    mp_process_tree(msghandler, tree->children, abdb, xbuf);
                    clean_xbuf(xbuf, { elem });
                    break;
                case AtomBondtypeXml::ATOMTYPE:
                    {
                        std::vector<AtomBondtypeXml> clean = {
                            AtomBondtypeXml::NAME, AtomBondtypeXml::INDEX
                        };
                        if (xmlFound(xbuf, clean))
                        {
                            abdb->back().atomtypes.push_back(xbuf->find(AtomBondtypeXml::NAME)->second);
                            clean_xbuf(xbuf, clean);
                        }
                        else
                        {
                            msghandler->msg(ACTStatus::Warning,
                                            gmx::formatString("Not all data present for determining an atomtype for %s", abdb->back().name.c_str()));
                        }
                    }
                    break;
                case AtomBondtypeXml::BONDTYPE:
                    {
                        std::vector<AtomBondtypeXml> clean = {
                            AtomBondtypeXml::AI, AtomBondtypeXml::AJ, AtomBondtypeXml::ORDER
                        };
                        if (xmlFound(xbuf, clean))
                        {
                            SimpleBond sb;
                            sb.ai    = xbuf_atoi(xbuf, AtomBondtypeXml::AI, true);
                            sb.aj    = xbuf_atoi(xbuf, AtomBondtypeXml::AJ, true);
                            sb.order = xbuf_atof(xbuf, AtomBondtypeXml::ORDER, true);
                            abdb->back().bonds.push_back(std::move(sb));
                            clean_xbuf(xbuf, clean);
                        }
                        else
                        {
                            msghandler->msg(ACTStatus::Warning,
                                            gmx::formatString("Not all data present for determining a bondtype for %s", abdb->back().name.c_str()));
                        }
                    }
                    break;
                default:
                    break;
                }
            }
        }
        tree = tree->next;
    }
}

void readAtomBondtypeDB(MsgHandler                     *msghandler,
                        const std::string              &dbname,
                        std::vector<AtomBondtypeEntry> *abdb)
{
    xmlDocPtr   doc;
    std::string mpfile;
    std::string database(dbname);

    if (database.empty())
    {
        database.assign("atom_bond.xml");
    }
    bool bAddCWD = true;
    bool bFatal  = false;
    auto mydb    = gmx::findLibraryFile(database, bAddCWD, bFatal);

    fillMaps();
    if (msghandler)
    {
        msghandler->msg(ACTStatus::Verbose, gmx::formatString("Opening %s\n", dbname.c_str()));
    }

    if ((doc = xmlParseFile(mydb.c_str())) == nullptr)
    {
        if (msghandler)
        {
            msghandler->msg(ACTStatus::Error,
                            gmx::formatString("Failed reading XML file %s. Run a syntax checker such as nsgmls.",
                                              mydb.c_str()));
        }
    }

    std::map<AtomBondtypeXml, std::string>  xbuf;
    mp_process_tree(msghandler, doc->children, abdb, &xbuf);
    if (msghandler)
    {
        msghandler->msg(ACTStatus::Verbose,
                        gmx::formatString("Read %zu special atom and bondtype entries from %s",
                                          abdb->size(), mydb.c_str()));
    }
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
