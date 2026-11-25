#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
"""!
@file
@brief Entries of the index of the manual that are made automatically:
the commands of alexandria, their command-line options and the variables
of CMake, wherever the text mentions them.

The commands and the options are read from the source code of ACT when the
manual is built, so that the index follows the program: the commands from
the calls of registerModule in src/act/alexandria/alex_modules.cpp, the
options from the strings "-name" in the code of src/act. An option that the
code does not know is not indexed, which keeps out words like non-bonded.
A boolean option written as -noname is indexed under -name. Options of one
letter, which other programs such as conda have too, count only in a text
that also has alexandria or one of its commands, and in a literal block only
on a command line of alexandria (with the lines after it that start with a
dash). A variable of CMake is a name with an underscore after -D.

Each name gets one entry for each section that mentions it, which points to
the first paragraph of the section that does. Entries of concepts are
written in the text with the directive index of Sphinx.
"""
import os
import re
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

## The source code of ACT
SOURCE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src", "act")

## The file that registers the commands of alexandria
MODULES_FILE = os.path.join("alexandria", "alex_modules.cpp")

## Commands of the command-line framework, which are not modules
BUILTIN_COMMANDS = ("help",)

## Options of the command-line framework of every command
BUILTIN_OPTIONS = ("h", "version", "quiet", "nice")

## Commands that are also English words: only indexed after the word
## alexandria, or where they are the whole text of a literal, an emphasis
## or a paragraph (the table of commands)
AMBIGUOUS_COMMANDS = ("simulate", "minimize", "help", "b2")

## Languages of literal blocks in which the options of the program are
## written, as on the command line
SHELL_LANGUAGES = ("", "default", "none", "text", "bash", "sh", "shell", "console")

## Languages of literal blocks in which no command or option is indexed
SKIPPED_LANGUAGES = ("c++", "cpp", "c", "json", "xml")

## What the index says each kind of entry is
KIND_TEXT = {
    "command": "alexandria command",
    "option": "command-line option",
    "cmake": "CMake variable",
}


def read_commands(source_dir: str = SOURCE_DIR) -> List[str]:
    """!
    @brief The commands of alexandria.
    @param source_dir The directory src/act of the source code.
    @return The names, sorted.
    """
    with open(os.path.join(source_dir, MODULES_FILE), encoding="utf-8") as f:
        text = f.read()
    names = set(re.findall(r'registerModule\s*\(\s*\w+\s*,\s*[^,]+,\s*"([a-z][a-z0-9_]*)"', text))
    names.update(BUILTIN_COMMANDS)
    return sorted(names)


def read_options(source_dir: str = SOURCE_DIR) -> List[str]:
    """!
    @brief The command-line options of the programs of ACT.
    @param source_dir The directory src/act of the source code.
    @return The names without the dash, sorted.
    """
    names = set(BUILTIN_OPTIONS)
    pattern = re.compile(r'"-([a-zA-Z][a-zA-Z0-9_]*)"')
    for root, dirs, files in os.walk(source_dir):
        dirs[:] = [d for d in dirs if d != "tests"]
        for name in files:
            if name.endswith(".cpp"):
                with open(os.path.join(root, name), encoding="utf-8", errors="replace") as f:
                    names.update(pattern.findall(f.read()))
    return sorted(names)


class Terms:
    """!
    @brief Finds the commands, options and variables of CMake in a text.
    """

    def __init__(self, commands: Iterable[str], options: Iterable[str]) -> None:
        """!
        @brief The names to look for.
        @param commands The commands of alexandria.
        @param options The options, without the dash.
        """
        self.commands = set(commands)
        self.options = set(options)
        clear = sorted((c for c in self.commands if c not in AMBIGUOUS_COMMANDS), key=len, reverse=True)
        self.clear_commands = re.compile(r"(?<![\w.-])(" + "|".join(map(re.escape, clear)) + r")(?![\w-])") \
            if clear else None
        self.after_program = re.compile(r"\balexandria\s+([a-z][a-z0-9_]*)\b")
        self.option = re.compile(r"(?<![\w.\-/=+*])-([a-zA-Z][a-zA-Z0-9_]*)(?![\w-])")
        self.cmake = re.compile(r"(?<![\w-])-D([A-Za-z][A-Za-z0-9_]*)=")

    def option_name(self, word: str) -> Optional[str]:
        """!
        @brief The option a word after a dash is.
        @param word The word, such as shell_keep or noshell_keep.
        @return The name of the option, or None if it is none.
        """
        if word in self.options:
            return word
        if word.startswith("no") and word[2:] in self.options:
            return word[2:]
        return None

    def find(self, text: str, whole: bool = False, options: bool = True,
             short: Optional[bool] = None, cmake: bool = True) -> List[Tuple[str, str]]:
        """!
        @brief The names in a text.
        @param text The text.
        @param whole Whether the text is a literal, an emphasis or a paragraph
               of its own, in which a command that is an English word counts.
        @param options Whether to look for options of the program.
        @param short Whether options of one letter count; by default only if
               the text has alexandria or a command.
        @param cmake Whether to look for variables of CMake.
        @return The kinds and names, in the order of the text, each once.
        """
        found: List[Tuple[int, str, str]] = []
        for m in self.after_program.finditer(text):
            if m.group(1) in self.commands:
                found.append((m.start(1), "command", m.group(1)))
        if self.clear_commands:
            for m in self.clear_commands.finditer(text):
                found.append((m.start(1), "command", m.group(1)))
        if whole:
            first = text.strip().split()[:2]
            if first and first[0] in self.commands and (len(first) == 1 or first[1].startswith("-")):
                found.append((text.find(first[0]), "command", first[0]))
        if short is None:
            short = bool(found) or re.search(r"\balexandria\b", text) is not None
        if cmake:
            for m in self.cmake.finditer(text):
                if "_" in m.group(1):
                    found.append((m.start(1), "cmake", m.group(1)))
        if options:
            for m in self.option.finditer(text):
                name = self.option_name(m.group(1))
                if name and not m.group(1).startswith("D") and (short or len(name) > 1):
                    found.append((m.start(1), "option", name))
        result: List[Tuple[str, str]] = []
        for _, kind, name in sorted(found):
            if (kind, name) not in result:
                result.append((kind, name))
        return result


def command_lines(block: str) -> List[Tuple[str, bool]]:
    """!
    @brief The command lines of a literal block.
    @param block The text of the block.
    @return The lines, each with the lines after it that start with a dash,
            and whether it runs alexandria.
    """
    result: List[Tuple[str, bool]] = []
    for line in block.split("\n"):
        if result and line.lstrip().startswith("-"):
            result[-1] = (result[-1][0] + "\n" + line, result[-1][1])
        else:
            result.append((line, re.search(r"\balexandria\b", line) is not None))
    return result


def entry(kind: str, name: str) -> Tuple[str, str]:
    """!
    @brief The entry of the index for a name.

    The options come first in the index, under Symbols, as in the HTML so in
    the PDF, whose index makeindex sorts.
    @param kind command, option or cmake.
    @param name The name.
    @return The type of the entry and its text.
    """
    if kind == "command":
        return "pair", KIND_TEXT[kind] + "; " + name
    if kind == "option":
        return "single", "-" + name + " (" + KIND_TEXT[kind] + ")"
    return "single", name + " (" + KIND_TEXT[kind] + ")"


def setup(app: Any) -> dict:
    """!
    @brief The extension of Sphinx: a transform that adds the entries.
    @param app The Sphinx application.
    @return The metadata of the extension.
    """
    from docutils import nodes  # pylint: disable=import-outside-toplevel
    from sphinx import addnodes  # pylint: disable=import-outside-toplevel
    from sphinx.transforms import SphinxTransform  # pylint: disable=import-outside-toplevel

    terms = Terms(read_commands(), read_options())

    ## Nodes whose text is not searched
    skipped = (nodes.math, nodes.math_block, nodes.raw, nodes.comment, nodes.reference,
               nodes.citation_reference, nodes.footnote_reference, nodes.target,
               nodes.substitution_definition, addnodes.pending_xref, addnodes.toctree,
               addnodes.desc_signature, addnodes.index)
    ## Inline nodes whose whole text counts as a word of its own
    whole_inlines = (nodes.literal, nodes.emphasis, nodes.strong)

    def is_skipped(node: Any) -> bool:
        """!
        @brief Whether a text is in a node that is not searched.
        """
        parent = node.parent
        while parent is not None:
            if isinstance(parent, skipped):
                return True
            if isinstance(parent, nodes.literal_block) and parent.get("language", "") in SKIPPED_LANGUAGES:
                return True
            parent = parent.parent
        return False

    def anchor_of(node: Any) -> Tuple[Any, bool]:
        """!
        @brief The element before which the entry is put.
        @return The element, and whether it is a title, after which the
                entry is put instead.
        """
        current = node.parent
        while current is not None:
            if isinstance(current, nodes.title) and isinstance(current.parent, nodes.section):
                return current, True
            if isinstance(current, nodes.Body) and current.parent is not None \
               and not isinstance(current.parent, nodes.TextElement):
                return current, False
            current = current.parent
        return None, False

    def section_of(node: Any) -> str:
        """!
        @brief The identifier of the section a node is in.
        """
        current = node.parent
        while current is not None:
            if isinstance(current, nodes.section) and current["ids"]:
                return current["ids"][0]
            current = current.parent
        return ""

    class ActIndexTransform(SphinxTransform):
        """!
        @brief Adds the entries, before Sphinx collects those of the index.
        """

        ## Before SphinxDomains (850), which collects the entries
        default_priority = 840

        def apply(self, **kwargs: Any) -> None:
            """!
            @brief Find the names and put the entries in the document.
            """
            seen: Set[Tuple[str, str, str]] = set()
            anchors: Dict[int, Tuple[Any, bool, List[Tuple[str, str]]]] = {}
            for text in list(self.document.findall(nodes.Text)):
                if is_skipped(text):
                    continue
                parent = text.parent
                block = parent if isinstance(parent, nodes.literal_block) else None
                language = block.get("language", "") if block is not None else ""
                shell = block is None or language in SHELL_LANGUAGES
                whole = isinstance(parent, whole_inlines) or \
                    (isinstance(parent, nodes.paragraph) and len(parent.children) == 1)
                if block is not None and shell:
                    found = []
                    for line, program in command_lines(str(text)):
                        for name in terms.find(line, options=program, short=program):
                            if name not in found:
                                found.append(name)
                else:
                    found = terms.find(str(text), whole=whole, options=shell, cmake=shell)
                if not found:
                    continue
                section = section_of(text)
                anchor, after = anchor_of(text)
                if anchor is None:
                    continue
                for kind, name in found:
                    if (section, kind, name) in seen:
                        continue
                    seen.add((section, kind, name))
                    anchors.setdefault(id(anchor), (anchor, after, []))[2].append((kind, name))
            for anchor, after, names in anchors.values():
                target_id = "index-%s" % self.env.new_serialno("index")
                target = nodes.target("", "", ids=[target_id])
                index = addnodes.index(entries=[], inline=False)
                for kind, name in names:
                    entry_type, value = entry(kind, name)
                    index["entries"].append((entry_type, value, target_id, "", None))
                parent = anchor.parent
                position = parent.index(anchor) + (1 if after else 0)
                parent.insert(position, [target, index])

    app.add_transform(ActIndexTransform)
    return {"version": "1.0", "parallel_read_safe": True, "parallel_write_safe": True}
