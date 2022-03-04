/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "act_checksum.h"

#include <unistd.h>

#include <cstdio>
#include <cstdlib>

#include <string>

#include "gromacs/fileio/md5.h"
#include "gromacs/utility/stringutil.h"

#include "poldata.h"
#include "poldata_xml.h"

namespace alexandria
{

std::string computeCheckSum(const std::string &filename)
{
    // Read the file as binary!    
#define BUFFERSIZE 1000000
    unsigned char buffer[BUFFERSIZE];
    FILE *filp = fopen(filename.c_str(), "rb"); 
    // Now compute the md5 checksum
    md5_state_t pms;
    int bytes_read;
    gmx_md5_init(&pms);
    do
    {
        bytes_read = fread(buffer, sizeof(unsigned char), BUFFERSIZE, filp);
        if (bytes_read > 0)
        {
            gmx_md5_append(&pms, buffer, bytes_read);
        }
    }
    while (bytes_read > 0);
    fclose(filp);
    
    auto mysum = gmx_md5_finish(&pms);
    // Now convert the checksum into a hex string
    std::string newVersion;
    for(const auto &m : mysum)
    {
        newVersion += gmx::formatString("%02x", m);
    }
    return newVersion;
}

std::string poldataCheckSum(Poldata *pd)
{
    // Save the old checkSum
    std::string oldCheckSum = pd->checkSum();
    pd->setCheckSum("");
    // Save the old timeStamp
    std::string timeStamp = pd->timeStamp();
    pd->setTimeStamp("");
    // Write the document to a tmp file
    char tmpPoldata[] = "poldataTmpXXXXXX";
    int fileDescriptor = mkstemp(tmpPoldata);
    if (fileDescriptor >= 0)
    {
        close(fileDescriptor);
        writePoldata(tmpPoldata, pd, false);
        // Now compute the md5 checksum
        std::string checksum = computeCheckSum(tmpPoldata);
        // Restore old chekcsum and time stamp
        pd->setCheckSum(oldCheckSum);
        pd->setTimeStamp(timeStamp);
        // Delete the tmp file
        int errcode = std::remove(tmpPoldata);
        if (errcode != 0)
        {
            std::perror("Could not delete temporary file");
        }
        return checksum;
    }
    else
    {
        fprintf(stderr, "Could not create a temporary file for generating a checksum.\n");
        return std::string("");
    }
}

} // namespace alexandria
