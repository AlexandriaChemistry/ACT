#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, re #, csv
from collections import namedtuple
#from enum import Enum
from elements import *

def parse_formula(formula):
  form = dict()
  if type(formula) is str:
    if not AtomNameToAtomNumber("H") == 1:
      sys.exit("Something fishy with reading elements")

    Element  = namedtuple('Element', ['name', 'counter'])
    elements = [Element(*reg) for reg in re.findall(r'([A-Z][a-z]*)(\d*)', formula)]
    for element in elements:
      if StringIsElement(element.name):
        if not element.name in form:
          form[element.name] = 0
        if (len(element.counter) > 0):
          form[element.name] += int(element.counter)
        else:
          form[element.name] += 1
        # if element.counter else True
      else:
        print("%s has undefined element %s" % (formula, element))
        return {}
  else:
    sys.exit("The input is not a string")
  form = {k:int(v) for k, v in form.items()}
  return form

def mol_size(formula):
    elems = parse_formula(formula)
    return sum(elems.values())

def mol_weight(formula):
    elems = parse_formula(formula)
    mw    = 0.0
    supp  = { "C": 12.0107, "H": 1.0079, "N": 14.0067, "O": 15.9994, "F": 18.9984, "P": 30.9738, "S": 32.065, "Cl": 35.453, "Br": 79.904, "I": 126.904 }
    for ee in elems.keys():
        if ee in supp:
            mw += elems[ee]*supp[ee]
        else:
            mw += 1000
    return mw
        
def supported(formula):
  supp     = [ "C", "H", "N", "O", "F", "P", "S", "Cl", "Br", "I" ]
  elements = parse_formula(formula)
  for elem in elements.keys():
      if elem not in supp:
          return False
  return True

def metal(formula):
  metals = [ "Li", "Na", "K", "Rb", "Cs", "Fr", "Be", "Mg", "Ca", "Sr", "Ba",
             "Ra", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
             "Al", "Ga", "Ge", "Pb", "Tl", "Hg", "Te", "Cd", "In", "Y", "Zr",
             "Ag", "La", "Sb", "Bi", "Po" ]

  elements = parse_formula(formula)
  for elem in elements.keys():
      if elem in metals:
          return True
  return False

def cmp_form_low(form_a, form_b):
    ffa = form_a["form"]
    ffb = form_b["form"]

    if ffa and len(ffa) > 0 and ffb and len(ffb) > 0:
        fa = parse_formula(ffa)
        fb = parse_formula(ffb)

        if not fa or not fb:
            return 0

        diff = 0
        for check in [ "C", "H", "I", "Br", "Cl", "S", "P", "F", "O", "N" ]:
            if check in fa:
                if check in fb:
                    diff = fa[check] - fb[check]
                else:
                    diff = 1
            elif check in fb:
                diff = -1
            if diff != 0:
                return diff
        # This assumes formulae are correct and sorted from high to low
        # atomic number, but this is not enforced in parse_formula
        for elem in fa:
            if elem in fb:
                diff = fa[elem] - fb[elem]
            if diff != 0:
                return diff

    return 0

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

def sort_formula(mylist):
    global elements
    elements = get_elements()
    return sorted(mylist, key=cmp_to_key(cmp_form_low))
