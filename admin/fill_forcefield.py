#!/usr/bin/env python3
#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/AlexandriaChemistry
#
# Copyright (c) 2026, David van der Spoel, Marie-Madeleine Walz, Paul J. van Maaren
#
"""Fill the bonded gaps of an ACT force field from other force fields.

Only parameters that the target lacks are added, so molecules that the
target already describes keep their energies:
- bonds, angles and linear angles get their geometry (bond length, angle,
  position on a linear angle) from the first donor that has them, and
  the force constants that most entries of the target have,
- a linear angle that the target has in one direction only gets the
  other direction, with 1-a,
- proper dihedrals get zeros for all Fourier coefficients, the value that
  most dihedrals of the target have, and ntrain 0,
- each new bond gets a bond charge correction with the default values.
Angles and dihedrals that no donor has can be given with --angle and
--dihedral. Only bonded types whose atoms all have a bond type in the
target are added, and only from donors with the same potential. The new
entries are inserted as text at the end of their interaction, after which
alexandria edit_ff can write the force field in its own layout.

ACS-pg1 in src/act/testdata/forcefields was made with, in that directory,
    fill_forcefield.py ACS-pg1.xml filled.xml ACS-pg-vs2.xml ACM-pg1.xml \
        ACS-pg2.xml ACM-pg2.xml lj-pc.xml PC+GV-elec.xml ACS-g1.xml ACS-g4.xml \
        --angle c2_b~o2_b~h_b=107.0 \
        --dihedral c2_b:c2_b~o2_b~h_b --dihedral h_b~c2_b~o2_b~h_b \
        --dihedral h_b~o2_b~c2_b=o2_b --dihedral c3_b~c2_b~o2_b~h_b
    alexandria edit_ff -ff filled.xml -o ACS-pg1.xml -write
The C-O-H angle of acids and phenols is that of the QM structure of acetic
acid (src/act/testdata/molecules/acetic-acid-3-oep.log.pdb). The torsions are not trained,
so the force field gives poor torsional modes.
"""
import argparse
import collections
import re
import sys
import xml.etree.ElementTree as ET

SEP = re.compile(r"([~=#:])")


def split(ident):
    parts = SEP.split(ident)
    return parts[0::2], parts[1::2]


def reverse(ident):
    atoms, bonds = split(ident)
    out = atoms[-1]
    for a, b in zip(reversed(atoms[:-1]), reversed(bonds)):
        out += b + a
    return out


def interactions(root):
    return {it.get("type"): it for it in root.iter("interaction")}


def lists(it):
    return {pl.get("identifier"): {p.get("type"): p.attrib for p in pl.iter("parameter")}
            for pl in it.iter("parameterlist")}


def known(ident, have, canswap):
    if ident in have:
        return True
    return canswap and reverse(ident) in have


def param(attrs):
    keys = ["type", "unit", "value", "uncertainty", "minimum", "maximum",
            "ntrain", "mutability", "nonnegative"]
    return "      <parameter %s/>\n" % " ".join('%s="%s"' % (k, attrs[k]) for k in keys)


def default(it, ptype):
    """The most common attributes of a force constant in the target."""
    count = collections.Counter()
    for pl in it.iter("parameterlist"):
        for p in pl.iter("parameter"):
            if p.get("type") == ptype:
                count[tuple(sorted(p.attrib.items()))] += 1
    return dict(count.most_common(1)[0][0])


def fmt(x):
    return ("%.4f" % x).rstrip("0").rstrip(".")


def scaled(attrs, value, ntrain=None, uncertainty="0"):
    out = dict(attrs)
    out["value"] = fmt(value)
    out["minimum"] = fmt(0.5*value)
    out["maximum"] = fmt(2*value)
    out["uncertainty"] = uncertainty
    if ntrain is not None:
        out["ntrain"] = ntrain
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="Force field to fill")
    ap.add_argument("output", help="Filled force field")
    ap.add_argument("donors", nargs="+", help="Force fields to take parameters from, in order")
    ap.add_argument("--angle", action="append", default=[], metavar="ID=DEGREES",
                    help="Add an angle that no donor has, e.g. c2_b~o2_b~h_b=107.0")
    ap.add_argument("--dihedral", action="append", default=[], metavar="ID",
                    help="Add a proper dihedral that no donor has, with zeros")
    args = ap.parse_args()

    text = open(args.target, encoding="ISO-8859-1").read()
    root = ET.fromstring(text.encode("ISO-8859-1"))
    tint = interactions(root)
    btypes = {o.get("value") for pt in root.iter("particletype")
              for o in pt.iter("option") if o.get("key") == "bondtype"}
    donors = [ET.parse(d).getroot() for d in args.donors]

    extra = {"ANGLES": {}, "PROPER_DIHEDRALS": {}}
    for a in args.angle:
        ident, value = a.split("=")
        extra["ANGLES"][ident] = {"angle": {"value": value, "uncertainty": "0", "ntrain": "1"}}
    for ident in args.dihedral:
        extra["PROPER_DIHEDRALS"][ident] = {}

    def candidates(itype, function):
        """Identifiers of itype in the donors, first donor first, then the extra ones."""
        seen = {}
        for d in donors:
            dit = interactions(d).get(itype)
            if dit is None or dit.get("function") != function:
                continue
            for ident, params in lists(dit).items():
                atoms, _ = split(ident)
                if all(a in btypes for a in atoms) and ident not in seen \
                   and reverse(ident) not in seen:
                    seen[ident] = params
        for ident, params in extra.get(itype, {}).items():
            if ident not in seen and reverse(ident) not in seen:
                seen[ident] = params
        return seen

    added = {}
    bondlength = {}
    for ident, params in lists(tint["BONDS"]).items():
        bondlength[ident] = float(params["bondlength"]["value"])

    def blength(a, sep, b):
        for ident in (a + sep + b, b + sep + a):
            if ident in bondlength:
                return bondlength[ident]
        return None

    for itype in ["BONDS", "ANGLES", "LINEAR_ANGLES", "PROPER_DIHEDRALS"]:
        it = tint[itype]
        canswap = it.get("canswap") != "false"
        have = lists(it)
        function = it.get("function")
        new = []
        for ident, dparams in candidates(itype, function).items():
            if known(ident, have, canswap):
                continue
            # A linear angle can't be swapped, but its reverse with 1-a is the same.
            if itype == "LINEAR_ANGLES" and reverse(ident) in have:
                continue
            body = ""
            if itype == "BONDS":
                if "bondlength" not in dparams:
                    continue
                b = dparams["bondlength"]
                for p in ["D0", "De", "beta"]:
                    body += param(default(it, p))
                body += param(scaled(default(it, "bondlength"), float(b["value"]),
                                     b["ntrain"], b["uncertainty"]))
                bondlength[ident] = float(b["value"])
            elif itype == "ANGLES":
                a = dparams.get("angle")
                if a is None:
                    continue
                angle = dict(default(it, "angle"))
                for k in ["value", "uncertainty", "ntrain"]:
                    angle[k] = a[k]
                v = float(a["value"])
                angle["minimum"] = fmt(0.5*v)
                angle["maximum"] = "180"
                body += param(angle) + param(default(it, "kt"))
            elif itype == "LINEAR_ANGLES":
                a = dparams.get("a")
                atoms, bonds = split(ident)
                b1 = blength(atoms[0], bonds[0], atoms[1])
                b2 = blength(atoms[1], bonds[1], atoms[2])
                if a is None or b1 is None or b2 is None:
                    continue
                body += param(dict(a))
                body += param(default(it, "klin")) + param(default(it, "kublin"))
                body += param(scaled(default(it, "r13lin"), b1 + b2))
            else:
                for c in ["c0", "c1", "c2", "c3", "c4", "c5"]:
                    p = dict(default(it, c))
                    p["value"] = "0"
                    p["ntrain"] = "0"
                    body += param(p)
            new.append('    <parameterlist identifier="%s">\n%s    </parameterlist>\n' % (ident, body))
            have[ident] = None
        if itype == "LINEAR_ANGLES":
            # A linear angle can't be swapped, so the target lists both
            # directions. Add a missing one from the other, with 1-a.
            old = lists(it)
            for ident, params in old.items():
                rident = reverse(ident)
                if rident in have:
                    continue
                a = dict(params["a"])
                a["value"] = fmt(1 - float(a["value"]))
                a["minimum"], a["maximum"] = (fmt(1 - float(params["a"]["maximum"])),
                                              fmt(1 - float(params["a"]["minimum"])))
                body = param(a) + "".join(param(params[p]) for p in ["klin", "kublin", "r13lin"])
                new.append('    <parameterlist identifier="%s">\n%s    </parameterlist>\n' % (rident, body))
                have[rident] = None
        added[itype] = new

    # The charge model needs a bond correction for each bond, between zeta types.
    zeta = {}
    for pt in root.iter("particletype"):
        opts = {o.get("key"): o.get("value") for o in pt.iter("option")}
        if "bondtype" in opts and "zetatype" in opts:
            zeta.setdefault(opts["bondtype"], opts["zetatype"])
    if "BONDCORRECTIONS" in tint:
        it = tint["BONDCORRECTIONS"]
        have = lists(it)
        new = []
        for ident in list(lists(tint["BONDS"])) + [re.search(r'identifier="([^"]*)"', b).group(1)
                                                    for b in added["BONDS"]]:
            atoms, bonds = split(ident)
            if atoms[0] not in zeta or atoms[1] not in zeta:
                continue
            za, zb = sorted([zeta[atoms[0]], zeta[atoms[1]]])
            bcc = za + bonds[0] + zb
            if bcc in have or reverse(bcc) in have:
                continue
            chi = dict(default(it, "delta_chi"))
            if za == zb:
                chi.update(value="0", minimum="0", maximum="0", mutability="Fixed", nonnegative="yes")
            else:
                chi.update(value="0", minimum="-8", maximum="8", mutability="Bounded", nonnegative="no")
            body = param(chi) + param(default(it, "delta_eta"))
            new.append('    <parameterlist identifier="%s">\n%s    </parameterlist>\n' % (bcc, body))
            have[bcc] = None
        added["BONDCORRECTIONS"] = new

    # Insert the new parameter lists at the end of each interaction.
    for itype, new in added.items():
        if not new:
            continue
        start = re.search(r'<interaction type="%s"[^>]*>' % itype, text)
        end = text.index("  </interaction>", start.end())
        text = text[:end] + "".join(new) + text[end:]
        print("%-18s %5d added" % (itype, len(new)), file=sys.stderr)

    open(args.output, "w", encoding="ISO-8859-1").write(text)


if __name__ == "__main__":
    main()
