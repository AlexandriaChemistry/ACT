#!/usr/bin/env python3
#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/AlexandriaChemistry
#
# Copyright (c) 2026, David van der Spoel, Marie-Madeleine Walz, Paul J. van Maaren
#
"""Clean up the ACT code with clang-tidy, the "tidy" target of make.

clang-tidy is run on the source files of src/act in the compilation
database of a build directory (compile_commands.json) with the checks of
admin/clang-tidy-fix.yaml. The fixes of all files are collected and
applied here, in the layout of the ACT code:
- const is written in front of the type (const auto x, not auto const x),
  and the padding of aligned declarations takes up the const where it can,
- braces added around statements are placed on lines of their own,
  indented by 4 spaces as the rest of the code,
- "= default" follows the parameter list on the same line.
Fixes that make scalar variables (built-in types, enumerations and
pointers) const are left out, unless --const-scalars is given.
Fixes that overlap with another fix are left for the next run, so run
make tidy until it reports no more fixes.

Only code that is compiled in the build directory is seen by clang-tidy,
so code for other configurations (e.g. SYCL or no MPI) is left as is.
Check the result with git diff and rebuild before committing.
"""

import argparse
import collections
import concurrent.futures
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

try:
    import yaml
except ImportError:
    sys.exit("clang_tidy_fix.py needs PyYAML (pip install pyyaml)")

# The built-in scalar types, as clang prints them
BUILTIN_SCALARS = {
    "bool", "char", "signed char", "unsigned char", "wchar_t", "char8_t", "char16_t", "char32_t",
    "short", "unsigned short", "int", "unsigned int", "unsigned", "long", "unsigned long",
    "long long", "unsigned long long", "__int128", "unsigned __int128",
    "float", "double", "long double", "std::nullptr_t", "nullptr_t",
}

CONST_MESSAGE = re.compile(r"variable '[^']*' of type '([^']*)'(?: \(aka '([^']*)'\))? can be declared 'const'")
ENUM_DECLARATION = re.compile(r"\benum\s+(?:class\s+|struct\s+)?([A-Za-z_]\w*)\s*(?::[^{;]*)?[{;]")
# Specifiers that stay in front of const
LEADING_SPECIFIERS = re.compile(rb"((?:static|thread_local|inline|constinit)\s+)*")
# The start of a line with a statement that may get braces
CONTROL_STATEMENT = re.compile(rb"\s*(\}\s*)?(else\s+if\b|if\b|for\b|while\b|else\b)")
# Fixes that are applied first when fixes overlap
PRIORITY = ["performance-unnecessary-copy-initialization", "performance-for-range-copy"]
# Width of "const "
CONST_WIDTH = 6


def enumerations(source_dir):
    """The names of all enumerations in the source code."""
    names = set()
    for top in ("src/act", "src/gromacs"):
        for root, _, files in os.walk(os.path.join(source_dir, top)):
            for name in files:
                if name.endswith((".h", ".cpp", ".hpp")):
                    with open(os.path.join(root, name), errors="replace") as handle:
                        names.update(ENUM_DECLARATION.findall(handle.read()))
    return names


def is_scalar(type_name, enums):
    """Whether a variable of this type (as clang prints it) is a scalar."""
    # A reference is judged by what it refers to
    name = type_name.strip().rstrip("&").strip()
    if name.endswith("*"):
        return True
    for qualifier in ("const ", "volatile "):
        name = name.replace(qualifier, "")
    name = name.strip()
    if name in BUILTIN_SCALARS:
        return True
    # Enumerations, also when written with their namespace or class
    return name.split("::")[-1] in enums


def is_scalar_const_fix(message, enums):
    """Whether a fix of misc-const-correctness concerns a scalar variable."""
    match = CONST_MESSAGE.search(message)
    return bool(match) and any(t and is_scalar(t, enums) for t in match.groups())


def line_start(text, offset):
    """The offset of the start of the line that contains offset."""
    return text.rfind(b"\n", 0, offset) + 1


def indentation(text, offset):
    """The indentation of the line that contains offset."""
    start = line_start(text, offset)
    end = start
    while end < len(text) and text[end:end + 1] in (b" ", b"\t"):
        end += 1
    return text[start:end]


def west_const(text, diag_offset, replacements):
    """Move the const that clang-tidy inserts before the variable name to the front
    of the declaration, keeping aligned declarations aligned where the padding
    allows it. None if not possible."""
    if len(replacements) != 1:
        return None
    (offset, length, new), = replacements
    if length != 0 or new.strip() != b"const":
        return None
    start = diag_offset + len(LEADING_SPECIFIERS.match(text, diag_offset).group(0))
    declaration = text[start:offset]
    if not declaration.strip() or b"\n" in declaration or b"(" in declaration or b"," in declaration:
        return None
    if declaration.lstrip().startswith(b"const "):
        return None
    fix = [(start, 0, b"const ")]
    # The padding between the type and the name of aligned declarations takes
    # up what it can of the const, keeping at least one space
    padding = len(declaration) - len(declaration.rstrip(b" "))
    removed = min(max(padding - 1, 0), CONST_WIDTH)
    if removed > 0:
        fix.append((offset - removed, removed, b""))
    shift = CONST_WIDTH - removed
    if shift > 0:
        # The rest of the declaration moves to the right, and so must the lines
        # that continue it where they are aligned with it: the lines in
        # parentheses or brackets, as arguments under an opening parenthesis,
        # and the lines outside brackets that start to the right of the = or (
        # after the name. Lines in braces, as the body of a lambda or the
        # elements of a list, are indented relative to the statement and stay
        # where they are.
        end = statement_end(text, offset)
        if end is None:
            return None
        after_name = re.compile(rb"[A-Za-z_]\w*\s*").match(text, offset)
        if after_name is None:
            return None
        aligned_column = after_name.end() - line_start(text, offset)
        for line in continuation_lines(text, offset, end):
            stack = continuation_lines.stack[line]
            if text[line:line + 1] in (b"\n", b"") or b"{" in stack:
                continue
            if stack or len(indentation(text, line)) > aligned_column:
                fix.append((line, 0, b" " * shift))
    return fix


def continuation_lines(text, offset, end):
    """The starts of the lines after the one of offset up to end, skipping
    strings and comments. For each, continuation_lines.stack holds the
    brackets that are open at its start."""
    starts = []
    continuation_lines.stack = {}
    stack = []
    i = offset
    while i < end:
        c = text[i:i + 1]
        if c in (b'"', b"'"):
            i += 1
            while i < end and text[i:i + 1] != c:
                i += 2 if text[i:i + 1] == b"\\" else 1
        elif text.startswith(b"//", i):
            i = text.find(b"\n", i) - 1
            if i < 0:
                break
        elif text.startswith(b"/*", i):
            i = text.find(b"*/", i) + 1
            if i <= 0:
                break
        elif c in (b"(", b"[", b"{"):
            stack.append(c)
        elif c in (b")", b"]", b"}"):
            if stack:
                stack.pop()
        elif c == b"\n":
            starts.append(i + 1)
            continuation_lines.stack[i + 1] = b"".join(stack)
        i += 1
    return starts


def statement_end(text, offset):
    """The offset of the semicolon that ends the statement in which offset lies,
    skipping brackets, strings and comments. None if not found."""
    depth = 0
    i = offset
    while i < len(text):
        c = text[i:i + 1]
        if c in (b'"', b"'"):
            i += 1
            while i < len(text) and text[i:i + 1] != c:
                i += 2 if text[i:i + 1] == b"\\" else 1
        elif text.startswith(b"//", i):
            i = text.find(b"\n", i)
            if i < 0:
                return None
        elif text.startswith(b"/*", i):
            i = text.find(b"*/", i)
            if i < 0:
                return None
            i += 1
        elif c in (b"(", b"[", b"{"):
            depth += 1
        elif c in (b")", b"]", b"}"):
            depth -= 1
            if depth < 0:
                return None
        elif c == b";" and depth == 0:
            return i
        i += 1
    return None


def braces(text, replacements):
    """Lay out the braces that clang-tidy adds around a statement on lines of their own.
    None if the fix does not have the expected form."""
    if len(replacements) != 2:
        return None
    (open_at, open_len, open_text), (close_at, close_len, close_text) = sorted(replacements)
    if open_len or close_len or open_text.strip() != b"{" or close_text.strip() != b"}":
        return None
    statement = text[open_at:close_at]
    if not statement.strip():
        return None
    # Only a statement that starts its line, as in "if (x) y;", not in one-line
    # functions such as "void f() { if (x) y; }", and that ends its line
    if not CONTROL_STATEMENT.match(text[line_start(text, open_at):open_at]):
        return None
    rest = text[close_at:text.find(b"\n", close_at) if b"\n" in text[close_at:] else len(text)].strip()
    if rest and not rest.startswith(b"//"):
        return None
    indent = indentation(text, open_at)
    lines = statement.strip().split(b"\n")
    body = [indent + b"    " + lines[0]]
    # Further lines of the statement keep their own indentation
    body.extend(lines[1:])
    new = b"\n" + indent + b"{\n" + b"\n".join(body) + b"\n" + indent + b"}"
    return [(open_at, close_at - open_at, new)]


def equals_default(text, replacements):
    """Put "= default;" on the line of the parameter list."""
    if len(replacements) != 1:
        return None
    (offset, length, new), = replacements
    if new.strip() != b"= default;":
        return None
    start = offset
    while start > 0 and text[start - 1:start] in (b" ", b"\t", b"\n"):
        start -= 1
    return [(start, offset + length - start, b" = default;")]


def realign(text, fix):
    """Keep the lines that continue the statement of a fix aligned with it:
    the arguments of a call or the parameters of a declaration move by what
    the fix changes before them on the first line, e.g. when push_back(Bond(
    becomes emplace_back( or virtual is removed."""
    first = min(offset for offset, _, _ in fix)
    first_line_end = text.find(b"\n", first)
    end = statement_end(text, first)
    if first_line_end < 0 or end is None or end < first_line_end:
        return fix
    delta = sum(len(new) - length for offset, length, new in fix if offset < first_line_end)
    first_column = first - line_start(text, first)
    if delta != 0:
        for line in continuation_lines(text, first, end):
            spaces = len(indentation(text, line))
            if (text[line:line + 1] in (b"\n", b"") or spaces <= first_column or
                b"{" in continuation_lines.stack[line]):
                continue
            if delta > 0:
                fix.append((line, 0, b" " * delta))
            else:
                fix.append((line, min(-delta, spaces - first_column - 1), b""))
    return fix


def emplace(text, replacements):
    """Remove the spaces inside the braces that modernize-use-emplace removes,
    such that push_back({ a, b }) becomes emplace_back(a, b)."""
    fix = []
    for offset, length, new in replacements:
        removed = text[offset:offset + length]
        if removed.endswith(b"{") and not new.endswith(b"{"):
            while text[offset + length:offset + length + 1] == b" ":
                length += 1
        elif new == b"" and removed == b"}":
            while offset > 0 and text[offset - 1:offset] == b" ":
                offset -= 1
                length += 1
        fix.append((offset, length, new))
    return realign(text, fix)


def reference(text, replacements):
    """Write the & that makes a variable a reference in front of the name, as
    in const auto &x, taking the place of a space of the padding of aligned
    declarations, and remove the line of a copy that is not used."""
    fix = []
    for offset, length, new in replacements:
        if new == b"" and length > 0:
            # A statement that is removed takes its line with it, if it is
            # the only statement on that line
            start = line_start(text, offset)
            removed = text[offset:offset + length]
            if text[start:offset].strip() == b"":
                if b"\n" in removed and removed[removed.rfind(b"\n"):].strip() == b"":
                    # The removal runs into the indentation of the next line
                    offset, length = start, offset + removed.rfind(b"\n") + 1 - start
                else:
                    end = text.find(b"\n", offset + length)
                    if end >= 0 and text[offset + length:end].strip() == b"":
                        offset, length = start, end + 1 - start
        elif length == 0 and new == b"&":
            spaces = 0
            while text[offset + spaces:offset + spaces + 1] == b" ":
                spaces += 1
            if spaces == 1:
                offset += 1
            elif spaces > 1:
                offset, length = offset + spaces - 1, 1
        fix.append((offset, length, new))
    return fix


def delete_null_pointer(text, replacements):
    """Replace the if statement around a delete, with its braces, by the
    delete statement on a line of its own, indented as the if was. None if
    the if does not have that form."""
    first = min(offset for offset, _, _ in replacements)
    last = max(offset + length for offset, length, _ in replacements)
    start = line_start(text, first)
    if not text[start:first].strip() == b"":
        return None
    match = re.compile(rb"\s*if\s*\([^;{}]*\)\s*\{?\s*(delete(?:\s*\[\s*\])?\s[^;{}]*;)\s*\}?").match(text, start)
    if match is None or match.end() < last:
        return None
    end = match.end()
    if text[end:end + 1] == b"\n":
        end += 1
    else:
        return None
    return [(start, end - start, indentation(text, start) + match.group(1) + b"\n")]


def layout(name, text, diag_offset, replacements):
    """The replacements of one fix in the layout of the ACT code, None to leave the fix out."""
    if name == "misc-const-correctness":
        return west_const(text, diag_offset, replacements)
    if name == "readability-braces-around-statements":
        return braces(text, replacements)
    if name == "modernize-use-equals-default":
        fix = equals_default(text, replacements)
        return fix if fix is not None else replacements
    if name == "modernize-use-emplace":
        return emplace(text, replacements)
    if name == "modernize-use-override":
        return realign(text, list(replacements))
    if name in ("performance-unnecessary-copy-initialization", "performance-for-range-copy"):
        return reference(text, replacements)
    if name == "readability-delete-null-pointer":
        return delete_null_pointer(text, replacements)
    return replacements


def overlaps(fix, taken):
    """Whether any replacement of fix touches a replacement that is already taken."""
    for offset, length, _ in fix:
        for start, end in taken:
            if offset == start or (offset < end and start < offset + length):
                return True
    return False


def selected_files(build_dir, source_dir, paths):
    """The source files of src/act in the compilation database, optionally limited to paths."""
    database = os.path.join(build_dir, "compile_commands.json")
    if not os.path.exists(database):
        sys.exit(f"No {database}, configure the build with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON")
    with open(database) as handle:
        entries = json.load(handle)
    act = os.path.join(source_dir, "src", "act") + os.sep
    files = set()
    for entry in entries:
        name = os.path.realpath(os.path.join(entry.get("directory", ""), entry["file"]))
        if name.startswith(act) and (not paths or any(name == p or name.startswith(p + os.sep) for p in paths)):
            files.add(name)
    return sorted(files)


def sysroot_args():
    """On macOS a clang-tidy that is not Apple's needs the SDK."""
    if sys.platform != "darwin":
        return []
    sdk = os.environ.get("SDKROOT")
    if not sdk:
        try:
            sdk = subprocess.run(["xcrun", "--show-sdk-path"], capture_output=True, text=True,
                                 check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return []
    return [f"--extra-arg=-isysroot{sdk}"] if sdk else []


def run_tidy(clang_tidy, build_dir, config, header_filter, extra, source, fixes):
    """Run clang-tidy on one file, the fixes are written to fixes."""
    command = [clang_tidy, "-p", build_dir, f"--config-file={config}", f"--header-filter={header_filter}",
               "--quiet", f"--export-fixes={fixes}"] + extra + [source]
    result = subprocess.run(command, capture_output=True, text=True)
    errors = [line for line in result.stdout.splitlines() + result.stderr.splitlines() if ": error:" in line]
    return source, errors


def collect_fixes(tmp, enums, const_scalars, allowed):
    """The fixes of all files: for each file a list of (check name, replacements),
    with the replacements of each fix as (offset, length, text), and the number of
    scalar const fixes that were left out."""
    contents = {}
    fixes = collections.defaultdict(list)
    seen = set()
    scalars = 0
    for name in sorted(os.listdir(tmp)):
        with open(os.path.join(tmp, name)) as handle:
            exported = yaml.safe_load(handle) or {}
        for diagnostic in exported.get("Diagnostics") or []:
            message = diagnostic.get("DiagnosticMessage", {})
            check = diagnostic.get("DiagnosticName", "")
            path = message.get("FilePath", "")
            replacements = message.get("Replacements") or []
            if not replacements or not allowed(path):
                continue
            if any(r["FilePath"] != path for r in replacements):
                continue
            if check == "misc-const-correctness" and not const_scalars and \
               is_scalar_const_fix(message.get("Message", ""), enums):
                scalars += 1
                continue
            key = (check, path, tuple((r["Offset"], r["Length"], r.get("ReplacementText", ""))
                                      for r in replacements))
            # A fix in a header is found once for each file that includes it
            if key in seen:
                continue
            seen.add(key)
            if path not in contents:
                with open(path, "rb") as handle:
                    contents[path] = handle.read()
            raw = [(r["Offset"], r["Length"], r.get("ReplacementText", "").encode()) for r in replacements]
            fix = layout(check, contents[path], message.get("FileOffset", 0), raw)
            fixes[path].append((check, fix))
    return contents, fixes, scalars


def apply_fixes(contents, fixes, dry_run):
    """Apply the fixes that do not overlap. Returns the number of fixes per check that were
    applied, left for the next run since they overlap, or left out since they cannot be
    written in the layout of the ACT code, and the number of changed files."""
    applied = collections.Counter()
    skipped = collections.Counter()
    unsupported = collections.Counter()
    changed = 0
    for path, file_fixes in sorted(fixes.items()):
        file_fixes.sort(key=lambda f: PRIORITY.index(f[0]) if f[0] in PRIORITY else len(PRIORITY))
        taken = []
        chosen = []
        for check, fix in file_fixes:
            if fix is None:
                unsupported[check] += 1
                continue
            if overlaps(fix, taken):
                skipped[check] += 1
                continue
            taken.extend((offset, offset + length) for offset, length, _ in fix)
            chosen.extend(fix)
            applied[check] += 1
        if chosen:
            changed += 1
        if chosen and not dry_run:
            text = contents[path]
            for offset, length, new in sorted(chosen, reverse=True):
                text = text[:offset] + new + text[offset + length:]
            with open(path, "wb") as handle:
                handle.write(text)
    return applied, skipped, unsupported, changed


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="*",
                        help="Files or directories (relative to the source directory) to limit the clean-up to, "
                        "default all of src/act; also from the environment variable ACT_TIDY_PATHS")
    parser.add_argument("--build-dir", default=".", help="Build directory with compile_commands.json")
    parser.add_argument("--source-dir", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        help="Top of the ACT source tree")
    parser.add_argument("--clang-tidy", default="clang-tidy", help="clang-tidy executable")
    parser.add_argument("--jobs", "-j", type=int, default=os.cpu_count(), help="Number of clang-tidy processes")
    parser.add_argument("--const-scalars", action="store_true",
                        help="Also make scalar variables (built-in types, enumerations, pointers) const")
    parser.add_argument("--dry-run", action="store_true", help="Only report the fixes, change nothing")
    # Options may also come from the environment, e.g. make tidy ACT_TIDY_ARGS=--dry-run
    args = parser.parse_args(sys.argv[1:] + shlex.split(os.environ.get("ACT_TIDY_ARGS", "")))

    source_dir = os.path.realpath(args.source_dir)
    build_dir = os.path.realpath(args.build_dir)
    if not shutil.which(args.clang_tidy):
        sys.exit(f"{args.clang_tidy} not found, set it with cmake -DACT_CLANG_TIDY=...")
    paths = [os.path.realpath(os.path.join(source_dir, p))
             for p in (args.paths or os.environ.get("ACT_TIDY_PATHS", "").split())]
    files = selected_files(build_dir, source_dir, paths)
    if not files:
        sys.exit("No source files of src/act selected")
    # Headers are changed only within the selected paths
    roots = paths or [os.path.join(source_dir, "src", "act")]
    header_filter = "|".join(re.escape(p) for p in roots)

    def allowed(path):
        path = os.path.realpath(path)
        return any(path == p or path.startswith(p + os.sep) for p in roots)

    config = os.path.join(source_dir, "admin", "clang-tidy-fix.yaml")
    extra = sysroot_args()
    print(f"Running clang-tidy on {len(files)} files with {args.jobs} processes")
    with tempfile.TemporaryDirectory(prefix="act-tidy-") as tmp:
        failed = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            jobs = [pool.submit(run_tidy, args.clang_tidy, build_dir, config, header_filter, extra, source,
                                os.path.join(tmp, f"tidy{i}.yaml"))
                    for i, source in enumerate(files)]
            for done, job in enumerate(concurrent.futures.as_completed(jobs), 1):
                source, errors = job.result()
                if errors:
                    failed[source] = errors
                if done % 50 == 0 or done == len(files):
                    print(f"  {done}/{len(files)} files")
        contents, fixes, scalars = collect_fixes(tmp, enumerations(source_dir), args.const_scalars, allowed)

    applied, skipped, unsupported, changed = apply_fixes(contents, fixes, args.dry_run)
    print("Fixes" + (" that would be applied:" if args.dry_run else " applied:"))
    for check, count in applied.most_common():
        print(f"  {count:6d} {check}")
    for check, count in skipped.most_common():
        print(f"  {count:6d} {check} left for the next run, since they overlap with other fixes")
    for check, count in unsupported.most_common():
        print(f"  {count:6d} {check} left out, since they cannot be written in the ACT layout")
    if scalars:
        print(f"  {scalars:6d} misc-const-correctness of scalar variables left out")
    for source, errors in sorted(failed.items()):
        print(f"Compiler errors in {os.path.relpath(source, source_dir)}, its fixes may be incomplete:")
        for line in errors[:3]:
            print(f"    {line}")
    if args.dry_run or not applied:
        print("Nothing changed")
    else:
        print(f"Changed {changed} files, check them with git diff and rebuild")
    return 0


if __name__ == "__main__":
    sys.exit(main())
