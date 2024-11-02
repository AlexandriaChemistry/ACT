#!/usr/bin/env python3

import os, re

def extract_todo(topdir:str)->list:
    todo_com = []
    comment = r'\btodo\b(.*)'
    comm2   = r'\bfixme\b(.*)'
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
                        lline = line.lower()
                        match = re.search(comment, lline)
                        if not match:
                            match = re.search(comm2, lline)
                        if match:
                            todo_com.append((filepath, i+1, match.group(1).strip()))

    return todo_com

with open("TODO.md", "w", encoding='utf-8') as outf:
    outf.write("## Open questions in the ACT code\n")
    outf.write("To update this list, run the get_todo.py script in the ACT root dir.\n")
    os.chdir("src/act")
    for ttt in extract_todo("."):
        outf.write("+ file: %s line: %d %s\n" % ( ttt[0], ttt[1], ttt[2] ))
