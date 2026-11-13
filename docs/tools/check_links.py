#!/usr/bin/env python3
"""Check relative links and #anchors in the Markdown documentation.

Usage (from the repository root):  python docs/tools/check_links.py [files or folders]

Anchors are computed the way GitHub and Codeberg/Forgejo do it: lower case,
punctuation removed, spaces turned into dashes, duplicates numbered -1, -2, ...
Explicit <a id="..."> / <a name="..."> anchors are also accepted.
Exit code 1 if a broken link is found.
"""
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def slug(text):
    text = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", text)          # images
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)       # links -> text
    text = re.sub(r"<[^>]+>", "", text)                         # html tags
    text = text.replace("`", "").strip().lower()
    text = re.sub(r"[^\w\- ]", "", text)
    return text.replace(" ", "-")


def anchors(path, cache={}):
    if path not in cache:
        found, seen, in_code = set(), {}, False
        for line in open(path, encoding="utf-8"):
            if line.lstrip().startswith("```"):
                in_code = not in_code
                continue
            if in_code:
                continue
            m = re.match(r"^(#{1,6})\s+(.*?)\s*#*\s*$", line)
            if m:
                s = slug(m.group(2))
                n = seen.get(s, 0)
                found.add(s if n == 0 else "%s-%d" % (s, n))
                seen[s] = n + 1
            for a in re.findall(r'<a\s+(?:id|name)="([^"]+)"', line):
                found.add(a)
        cache[path] = found
    return cache[path]


def links(path):
    in_code = False
    for no, line in enumerate(open(path, encoding="utf-8"), 1):
        if line.lstrip().startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        line = re.sub(r"`[^`]*`", "", line)
        for m in re.finditer(r"\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)", line):
            yield no, m.group(1)
        for m in re.finditer(r'(?:src|href)="([^"]+)"', line):
            yield no, m.group(1)


def main(targets):
    files = []
    for t in targets or [os.path.join(REPO, "README.md"), os.path.join(REPO, "docs")]:
        if os.path.isdir(t):
            for dp, _, fs in os.walk(t):
                files += [os.path.join(dp, f) for f in fs if f.endswith(".md")]
        else:
            files.append(t)
    bad = 0
    for f in sorted(files):
        for no, target in links(f):
            if re.match(r"^[a-z]+:", target) or target.startswith("//"):
                continue
            page, _, anchor = target.partition("#")
            dest = os.path.normpath(os.path.join(os.path.dirname(f), page)) if page else f
            if not os.path.exists(dest):
                print("%s:%d: missing file %s" % (os.path.relpath(f, REPO), no, target))
                bad += 1
            elif anchor and dest.endswith(".md") and anchor not in anchors(dest):
                print("%s:%d: missing anchor %s" % (os.path.relpath(f, REPO), no, target))
                bad += 1
    print("%d files checked, %d broken links" % (len(files), bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
