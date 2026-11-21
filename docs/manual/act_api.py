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
@brief The reference of the Python package act in the manual, made from
the documentation in the code when the manual is built.

The directive act-api of the manual imports the package and describes its
functions and classes with the directives of the Python domain of Sphinx.
The classes of ACT come from the compiled module act._core, whose
documentation is in python/src/core.cpp; the rest of the package documents
itself with docstrings in the style of Doxygen, which are converted here.

The package is taken from the directory in the environment variable
ACT_PYTHON_PATH, which the build of the manual sets to the package of the
build tree, or else from python/ of the source tree, which has no classes.
"""
import importlib
import inspect
import os
import re
import sys
from typing import Any, List, Optional, Tuple

## The directory of the package act in the source tree
SOURCE_PACKAGE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python")

## Dunder methods of the classes that are part of the interface
DOCUMENTED_DUNDERS = ("__getitem__", "__len__")


def escape(text: str) -> str:
    """!
    @brief Escape the characters that reStructuredText would take for markup.
    @param text Plain text.
    @return The text, which renders as it is.
    """
    for c in ("\\", "*", "`", "|", "_"):
        text = text.replace(c, "\\" + c)
    return text


def doxygen_to_rst(doc: str) -> List[str]:
    """!
    @brief Convert a docstring in the style of Doxygen to reStructuredText.

    @@brief starts the text, @@param, @@return and @@throws become fields,
    @@code{.py} ... @@endcode a code block, and @@p name a name in code.
    @param doc The docstring, which may start with an exclamation mark.
    @return The lines.
    """
    lines = inspect.cleandoc(doc or "").split("\n")
    if lines and lines[0].strip() == "!":
        lines = lines[1:]
    out: List[str] = []
    field: Optional[str] = None
    in_code = False
    for line in lines:
        stripped = line.strip()
        if in_code:
            if stripped.startswith("@endcode"):
                in_code = False
                out.append("")
            else:
                out.append("   " + line)
            continue
        if stripped.startswith("@code"):
            field = None
            out.extend(["", ".. code-block:: python", ""])
            in_code = True
            continue
        if stripped.startswith("@file"):
            continue
        tag = re.match(r"@(brief|param|return|returns|throws|throw)\b\s*(.*)", stripped)
        if tag:
            kind, rest = tag.group(1), tag.group(2)
            # A list of fields starts after an empty line
            if kind != "brief" and field is None and out and out[-1]:
                out.append("")
            if kind == "brief":
                field = None
                out.append(inline(rest))
            elif kind == "param":
                name, _, text = rest.partition(" ")
                field = ":param %s: " % name
                out.append(field + inline(text.strip()))
            elif kind in ("return", "returns"):
                field = ":returns: "
                out.append(field + inline(rest))
            else:
                name, _, text = rest.partition(" ")
                field = ":raises %s: " % name
                out.append(field + inline(text.strip()))
            continue
        if not stripped:
            field = None
            out.append("")
        elif field:
            # The continuation of a field
            out.append("   " + inline(stripped))
        else:
            out.append(inline(stripped))
    return out


def inline(text: str) -> str:
    """!
    @brief Escape a line of a docstring, with @@p name as a name in code.
    @param text The line.
    @return The line in reStructuredText.
    """
    parts = re.split(r"@p\s+(\S+)", text)
    out = ""
    for i, part in enumerate(parts):
        out += ("``%s``" % part) if i % 2 == 1 else escape(part)
    return out


def plain_to_rst(doc: str) -> List[str]:
    """!
    @brief Convert the plain documentation of the compiled module.
    @param doc The text.
    @return The lines.
    """
    return [inline(line) for line in inspect.cleandoc(doc or "").split("\n")]


def clean_signature(signature: str) -> str:
    """!
    @brief The signature of a function of the compiled module as the manual shows it.
    @param signature The signature that pybind11 writes, such as
    minimize(self: act._core.Molecule, *, max_iter: int = None) -> act._core.Minimization
    @return The signature without self and with act for act._core.
    """
    signature = re.sub(r"self: [\w.]+(, )?", "", signature)
    # The types that pybind11 3 writes for what Python passes as numbers,
    # lists and files
    for long, short in (("typing.SupportsFloat | typing.SupportsIndex", "float"),
                        ("typing.SupportsInt | typing.SupportsIndex", "int"),
                        ("collections.abc.Sequence", "Sequence"),
                        ("os.PathLike | str | bytes", "str | os.PathLike"),
                        ("act._core.", "act.")):
        signature = signature.replace(long, short)
    return signature


def compiled_doc(obj: Any, name: str) -> Tuple[Optional[str], List[str]]:
    """!
    @brief The signature and the documentation of a function of the compiled module.
    @param obj  The function.
    @param name Its name.
    @return The signature, None if the documentation has none, and the lines of the text.
    """
    doc = inspect.cleandoc(getattr(obj, "__doc__", None) or "")
    first, _, rest = doc.partition("\n")
    if first.startswith(name + "("):
        return clean_signature(first), plain_to_rst(rest)
    return None, plain_to_rst(doc)


def python_signature(obj: Any, name: str) -> str:
    """!
    @brief The signature of a function or method written in Python.
    @param obj  The function.
    @param name Its name.
    @return The signature without self.
    """
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        return name + "(...)"

    def annotation(value: Any) -> str:
        """The annotation as written, also when the annotations are strings."""
        return value if isinstance(value, str) else inspect.formatannotation(value)

    parts = []
    star = False
    for p in sig.parameters.values():
        if p.name in ("self", "cls"):
            continue
        text = p.name
        if p.kind == p.VAR_POSITIONAL:
            text, star = "*" + text, True
        elif p.kind == p.VAR_KEYWORD:
            text = "**" + text
        elif p.kind == p.KEYWORD_ONLY and not star:
            parts.append("*")
            star = True
        if p.annotation is not p.empty:
            text += ": " + annotation(p.annotation)
        if p.default is not p.empty:
            text += " = " + repr(p.default)
        parts.append(text)
    result = "%s(%s)" % (name, ", ".join(parts))
    if sig.return_annotation is not sig.empty:
        result += " -> " + annotation(sig.return_annotation)
    return result.replace("'", "")


def indented(lines: List[str], level: int = 1) -> List[str]:
    """!
    @brief Indent lines for the content of a directive.
    @param lines The lines.
    @param level The number of levels of three spaces.
    @return The indented lines.
    """
    pad = "   " * level
    return [(pad + line) if line else "" for line in lines]


def describe_function(name: str, obj: Any, directive: str = "function", level: int = 0) -> List[str]:
    """!
    @brief Describe a function or a method.
    @param name      Its name.
    @param obj       The function.
    @param directive The directive of the Python domain, function, method or staticmethod.
    @param level     The indentation of the directive.
    @return The lines.
    """
    if inspect.isbuiltin(obj) or type(obj).__name__ in ("builtin_function_or_method", "instancemethod"):
        signature, text = compiled_doc(obj, name)
        signature = signature or name + "(...)"
    else:
        signature, text = python_signature(obj, name), doxygen_to_rst(obj.__doc__)
    lines = [".. py:%s:: %s" % (directive, signature), ""] + indented(text) + [""]
    return indented(lines, level)


def describe_class(name: str, cls: type) -> List[str]:
    """!
    @brief Describe a class with its methods and properties.
    @param name The name under which the package has it.
    @param cls  The class.
    @return The lines.
    """
    compiled = cls.__module__ == "act._core"
    text = plain_to_rst(cls.__doc__) if compiled else doxygen_to_rst(cls.__doc__)
    if issubclass(cls, BaseException):
        bases = [b.__name__ for b in cls.__bases__ if b is not object]
        return [".. py:exception:: %s" % name, ""] + indented(
            text + ([""] + ["Derived from %s." % ", ".join(bases)] if bases else [])) + [""]
    lines = [".. py:class:: %s" % name, ""] + indented(text) + [""]
    for member, raw in sorted(vars(cls).items()):
        if member.startswith("_") and member not in DOCUMENTED_DUNDERS:
            continue
        value = getattr(cls, member)
        if isinstance(raw, property) or type(raw).__name__ == "property":
            doc = raw.__doc__ or ""
            text = plain_to_rst(doc) if compiled else doxygen_to_rst(doc)
            lines += indented([".. py:property:: %s" % member, ""] + indented(text) + [""])
        elif isinstance(raw, staticmethod) or (compiled and inspect.isbuiltin(value)
                                                and getattr(value, "__self__", None) is cls):
            lines += describe_function(member, value, "staticmethod", 1)
        elif callable(value):
            lines += describe_function(member, value, "method", 1)
    return lines


def load_package() -> Any:
    """!
    @brief Import the package act from where the build of the manual says.
    @return The package.
    """
    path = os.environ.get("ACT_PYTHON_PATH") or SOURCE_PACKAGE_DIR
    sys.path.insert(0, os.path.abspath(path))
    for name in [m for m in sys.modules if m == "act" or m.startswith("act.")]:
        del sys.modules[name]
    try:
        return importlib.import_module("act")
    finally:
        sys.path.pop(0)


def api_lines(act: Any) -> List[str]:
    """!
    @brief The reference of the package as reStructuredText.
    @param act The package.
    @return The lines.
    """
    lines = [".. py:module:: act", ""]
    functions, classes, errors, objects = [], [], [], []
    names = list(act.__all__)
    if act.classes_available():
        names += list(act._CLASSES)  # pylint: disable=protected-access
    for name in sorted(set(names)):
        obj = getattr(act, name)
        if inspect.isclass(obj):
            (errors if issubclass(obj, BaseException) else classes).append((name, obj))
        elif callable(obj):
            functions.append((name, obj))
        else:
            objects.append((name, obj))
    lines += ["Functions", "---------", ""]
    for name, obj in functions:
        lines += describe_function(name, obj)
    for name, obj in objects:
        lines += [".. py:data:: %s" % name, ""] + indented(doxygen_to_rst(type(obj).__doc__)) + [""]
    lines += ["Classes", "-------", ""]
    if not act.classes_available():
        lines += ["This manual was built without the classes of ACT, which come with",
                  "ACT built with -DACT_PYTHON=ON or with the conda package alexandria-python.", ""]
    for name, obj in classes:
        lines += describe_class(name, obj)
    lines += ["Errors", "------", ""]
    for name, obj in errors:
        lines += describe_class(name, obj)
    return lines


def setup(app: Any) -> dict:
    """!
    @brief The extension of Sphinx: the directive act-api.
    @param app The Sphinx application.
    @return The metadata of the extension.
    """
    from docutils.statemachine import StringList  # pylint: disable=import-outside-toplevel
    from sphinx.util.docutils import SphinxDirective  # pylint: disable=import-outside-toplevel
    from sphinx.util.nodes import nested_parse_with_titles  # pylint: disable=import-outside-toplevel
    from docutils import nodes  # pylint: disable=import-outside-toplevel
    from sphinx.errors import ExtensionError  # pylint: disable=import-outside-toplevel

    class ActApiDirective(SphinxDirective):
        """!
        @brief Insert the reference of the package act.
        """

        def run(self) -> list:
            """!
            @brief Make the nodes of the reference.
            @return The nodes.
            """
            act = load_package()
            # The build of the manual gives the package of the build tree
            # when ACT was built with the classes; that it has none then
            # means that Sphinx runs another Python than the module was
            # built for, which would leave the classes out unnoticed
            if os.environ.get("ACT_PYTHON_PATH") and not act.classes_available():
                raise ExtensionError("The package act in %s has no classes: Sphinx (%s) has to run the "
                                     "Python that the module act._core was built for" %
                                     (os.environ["ACT_PYTHON_PATH"], sys.executable))
            node = nodes.section()
            node.document = self.state.document
            nested_parse_with_titles(self.state, StringList(api_lines(act), "act-api"), node)
            return node.children

    app.add_directive("act-api", ActApiDirective)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
