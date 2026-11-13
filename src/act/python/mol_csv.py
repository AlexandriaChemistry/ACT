#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# Copyright (C) 2026
#
# Developers:
#             Mohammad Mehdi Ghahremanpour,
#             Paul J. van Maaren,
#             David van der Spoel (Project leader)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
#
"""The molecules of alexandria.csv as JSON in UTF-8.

alexandria.csv has one molecule per line, with the columns separated by
'|', in the layout of version 3 of mol_csv_api.Molecule.from_row. Opening
and saving it in a spreadsheet program may change its encoding, therefore
the data can be kept as JSON instead, always written as UTF-8 with the
characters themselves (for instance a plus-minus sign) rather than escapes.
"""

import json

# The columns of alexandria.csv, version 3, in their order
COLUMNS = [ "alexandria", "parentid", "iupac", "formula", "charge",
            "multiplicity", "weight", "cas", "csid", "pubchem", "inchi",
            "inchikey", "classification", "synonyms", "filename",
            "point_group", "symmetry_number", "rotatable_bonds",
            "html_name", "latex_name", "author", "timestamp" ]
# An optional last column
COMMENT  = "comment"
# Columns that are integers
INTEGERS = [ "alexandria", "parentid", "charge", "multiplicity",
             "symmetry_number", "rotatable_bonds" ]
# Columns that are numbers with decimals, absent if empty
FLOATS   = [ "weight" ]
# Columns that are lists separated by ';'
LISTS    = [ "classification", "synonyms" ]

def fix_double_encoding(text:str) -> str:
    """Repair text that was encoded as UTF-8 twice, such as 'Â±' for '±'.

    This happens when a program reads UTF-8 as Latin-1 and saves it as
    UTF-8 again, for instance a spreadsheet. The text is returned as it
    is if it cannot be such a text, or if the repair does not give UTF-8.
    """
    if not any(c in text for c in "ÂÃ"):
        return text
    try:
        return text.encode("latin-1").decode("utf-8")
    except (UnicodeEncodeError, UnicodeDecodeError):
        return text

def fix_molecule_encoding(mol:dict) -> int:
    """Repair the strings of a molecule with fix_double_encoding.

    Returns the number of strings that changed.
    """
    nfix = 0
    for name, value in mol.items():
        if isinstance(value, str):
            fixed = fix_double_encoding(value)
            nfix += 1 if fixed != value else 0
            mol[name] = fixed
        elif isinstance(value, list):
            fixed = [ fix_double_encoding(v) for v in value ]
            nfix += sum(1 for a, b in zip(fixed, value) if a != b)
            mol[name] = fixed
    return nfix

def row_to_molecule(row:list) -> dict:
    """Convert the columns of a line of alexandria.csv to a dictionary.

    Integers and numbers become numbers, the classification and the
    synonyms lists of strings, and the other columns strings. An empty
    weight is left out, as is a comment that is not there.
    """
    if len(row) < len(COLUMNS):
        raise ValueError("A molecule needs %d columns, this line has %d: %s" %
                         ( len(COLUMNS), len(row), "|".join(row) ))
    mol = {}
    for name, value in zip(COLUMNS, row):
        if name in INTEGERS:
            mol[name] = int(value)
        elif name in FLOATS:
            if len(value) > 0:
                mol[name] = float(value)
        elif name in LISTS:
            mol[name] = [ v for v in value.split(";") if len(v) > 0 ]
        else:
            mol[name] = value
    if len(row) > len(COLUMNS) and len(row[len(COLUMNS)]) > 0:
        mol[COMMENT] = row[len(COLUMNS)]
    return mol

def molecule_to_row(mol:dict) -> list:
    """Convert a molecule of row_to_molecule back to the columns of alexandria.csv.

    Numbers are written such that they read back as the same number, a
    missing weight as an empty column, and the lists separated by ';'.
    """
    row = []
    for name in COLUMNS:
        value = mol.get(name)
        if value is None:
            row.append("")
        elif name in LISTS:
            row.append(";".join(value))
        elif name in FLOATS:
            row.append(repr(float(value)))
        else:
            row.append(str(value))
    if COMMENT in mol and mol[COMMENT]:
        row.append(mol[COMMENT])
    return row

def read_mol_csv(filename:str, delimiter:str="|", comment:str="#") -> list:
    """Read the molecules of alexandria.csv, in the order of the file.

    The file is read as UTF-8. Empty lines and lines that start with the
    comment character are skipped.
    """
    mols = []
    with open(filename, "r", encoding="utf-8") as inf:
        for line in inf:
            line = line.rstrip("\r\n")
            if len(line.strip()) == 0 or line.startswith(comment):
                continue
            mols.append(row_to_molecule(line.split(delimiter)))
    return mols

def write_mol_json(mols:list, filename:str, indent:int=None):
    """Write the molecules to a JSON file in UTF-8.

    The characters are written as they are, not as escapes, so that the
    file can be read and compared as text. By default each molecule is on
    a line of its own, such that a change of a molecule is a change of one
    line, as in the csv file; with indent each member is on its own line.
    """
    with open(filename, "w", encoding="utf-8") as outf:
        if indent is not None:
            json.dump(mols, outf, ensure_ascii=False, indent=indent)
            outf.write("\n")
            return
        outf.write("[\n")
        for i, mol in enumerate(mols):
            outf.write(json.dumps(mol, ensure_ascii=False))
            outf.write(",\n" if i + 1 < len(mols) else "\n")
        outf.write("]\n")

def read_mol_json(filename:str) -> list:
    """Read the molecules of a JSON file written by write_mol_json."""
    with open(filename, "r", encoding="utf-8") as inf:
        return json.load(inf)
