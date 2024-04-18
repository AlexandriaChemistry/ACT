#!/usr/bin/env python3

import os, re

def extract_todo(topdir:str)->list:
    todo_com = []
    comment = r'\bTODO\b:(.*)'
    comm2   = r'\bFixme\b:(.*)'
    for root, dirs, files in os.walk(top = topdir):
        for cpp in files:
            filepath = os.path.join(root, cpp)
            myextensions = [ ".cpp", ".h", ".py", ".cmakein", ".txt" ]
            found = False
            for myext in myextensions:
                if filepath.find(myext) == len(filepath) - len(myext):
                    found = True
            if found:
                with open(filepath, 'r', encoding='utf-8') as file:
                    lines = file.readlines()
                for i, line in enumerate(lines):
                    match = re.search(comment, line)
                    if not match:
                        match = re.search(comm2, line)
                    if match:
                        todo_com.append((filepath, i+1, match.group(1).strip()))

    return todo_com

with open("TODO.md", "w", encoding='utf-8') as outf:
    outf.write("## Open questions in the ACT code\n")
    os.chdir("src/act")
    for ttt in extract_todo("."):
        outf.write("+ file: %s line: %d %s\n" % ( ttt[0], ttt[1], ttt[2] ))
