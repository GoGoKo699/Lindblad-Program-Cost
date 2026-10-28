"""Check repository navigation, Markdown structure, and milestone records.

These dependency-free checks cover source structure. Hosted equation rendering
still requires inspection in GitHub's rendered view.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import unicodedata
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
FIELDS = (
    "Central question.",
    "Understanding target.",
    "Dependencies and story role.",
    "Conceptual connection.",
    "Source pointers.",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _escaped(text: str, index: int) -> bool:
    start = index
    while start and text[start - 1] == "\\":
        start -= 1
    return (index - start) % 2 == 1


def _mask(text: str) -> str:
    """Keep offsets and line numbers intact when excluding literal content."""
    return "".join("\n" if char == "\n" else " " for char in text)


def _validate_math(text: str, name: str, line_number: int) -> None:
    """Guard compatibility defects observed in the hosted GitHub renderer."""
    _require(not re.search(r"\\operatorname\b", text), f"{name}:{line_number}: GitHub rejects \\operatorname; use the supported \\mathrm spelling")
    _require(not re.search(r"[<>]", text), f"{name}:{line_number}: use \\lt or \\gt instead of raw angle brackets in math")


def _source_views(source: str, name: str) -> tuple[str, str, Counter]:
    """Return prose-only and math-masked views, plus delimiter counts."""
    lines = source.splitlines(keepends=True)
    visible = []
    fence = None
    counts = Counter()
    for line_number, line in enumerate(lines, 1):
        marker = re.match(r"^ {0,3}(`{3,}|~{3,})(.*)$", line.rstrip("\r\n"))
        if fence is not None:
            visible.append(_mask(line))
            if marker and marker[1][0] == fence[0] and len(marker[1]) >= fence[1] and not marker[2].strip():
                fence = None
            elif fence[3] == "math":
                _validate_math(line, name, line_number)
            continue
        if marker:
            fence = (marker[1][0], len(marker[1]), line_number, marker[2].strip())
            counts["fenced_blocks"] += 1
            counts["math_blocks"] += marker[2].strip() == "math"
            visible.append(_mask(line))
        else:
            visible.append(line)
    _require(fence is None, f"{name}:{fence[2] if fence else 0}: unclosed fenced block")

    text = "".join(visible)
    prose = list(text)
    math_masked = list(text)
    index = 0
    while index < len(text):
        if _escaped(text, index):
            index += 1
            continue
        if text.startswith("$`", index):
            end = text.find("`$", index + 2)
            line_number = text.count("\n", 0, index) + 1
            _require(end >= 0, f"{name}:{line_number}: unclosed $` inline-math delimiter")
            _require("$`" not in text[index + 2:end], f"{name}:{line_number}: nested $` inline-math delimiter")
            _validate_math(text[index + 2:end], name, line_number)
            end += 2
            prose[index:end] = _mask(text[index:end])
            # Retain a nonblank cell marker: an all-math first table cell is
            # still a cell when the table omits its optional outer pipes.
            math_masked[index:end] = "M" + _mask(text[index + 1:end])
            counts["inline_math"] += 1
            index = end
        elif text.startswith("`$", index):
            line_number = text.count("\n", 0, index) + 1
            raise AssertionError(f"{name}:{line_number}: unmatched `$ inline-math delimiter")
        elif text[index] == "`":
            run = re.match(r"`+", text[index:])[0]
            closing = re.search(r"(?<!`)" + re.escape(run) + r"(?!`)", text[index + len(run):])
            if closing is None:
                index += len(run)  # An unmatched ordinary backtick is literal Markdown.
            else:
                end = index + len(run) + closing.end()
                prose[index:end] = _mask(text[index:end])
                index = end
        else:
            if text[index] == "\\" and index + 1 < len(text) and text[index + 1] in "()[]":
                line_number = text.count("\n", 0, index) + 1
                raise AssertionError(f"{name}:{line_number}: use GitHub math delimiters instead of {text[index:index + 2]}")
            index += 1
    return "".join(prose), "".join(math_masked), counts


def _slug(heading: str) -> str:
    """GitHub's ordinary heading slug, including its duplicate suffixes below."""
    heading = re.sub(r"<[^>]*>", "", heading)
    heading = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", heading)
    heading = heading.lower()
    heading = "".join(char for char in heading if char in "-_ " or not unicodedata.category(char).startswith(("P", "S")))
    return heading.replace(" ", "-")


def _headings(text: str) -> tuple[set[str], list[tuple[int, int, str]]]:
    anchors = set()
    headings = []
    for line_number, line in enumerate(text.splitlines(), 1):
        match = re.match(r"^ {0,3}(#{1,6})[ \t]+(.*?)(?:[ \t]+#+[ \t]*)?$", line)
        if match:
            value = match[2].strip()
            base = _slug(value)
            anchor = base
            suffix = 0
            while anchor in anchors:
                suffix += 1
                anchor = f"{base}-{suffix}"
            anchors.add(anchor)
            headings.append((line_number, len(match[1]), value))
    return anchors, headings


def _destinations(text: str):
    """Yield inline-link and reference-definition destinations with offsets."""
    for match in re.finditer(r"(?<!!)\]\(", text):
        if _escaped(text, match.start()):
            continue
        start = match.end()
        while start < len(text) and text[start].isspace():
            start += 1
        if start < len(text) and text[start] == "<":
            end = start + 1
            while end < len(text) and (text[end] != ">" or _escaped(text, end)):
                end += 1
            if end < len(text):
                yield text[start + 1:end], start
            continue
        end, depth = start, 0
        while end < len(text):
            char = text[end]
            if not _escaped(text, end):
                if char == "(":
                    depth += 1
                elif char == ")":
                    if depth == 0:
                        break
                    depth -= 1
                elif char.isspace() and depth == 0:
                    break
            end += 1
        if end > start:
            yield re.sub(r"\\([()])", r"\1", text[start:end]), start
    for match in re.finditer(r"^ {0,3}\[[^\]\n]+\]:[ \t]*(<[^>\n]+>|\S+)", text, re.MULTILINE):
        yield match[1].strip("<>"), match.start(1)


def _cells(line: str) -> list[str]:
    positions = [i for i, char in enumerate(line) if char == "|" and not _escaped(line, i)]
    if not positions:
        return []
    starts = [-1] + positions
    ends = positions + [len(line)]
    cells = [line[start + 1:end].strip() for start, end in zip(starts, ends)]
    if line.lstrip().startswith("|"):
        cells.pop(0)
    if cells and line.rstrip().endswith("|") and not _escaped(line.rstrip(), len(line.rstrip()) - 1):
        cells.pop()
    return cells


def _tables(text: str, name: str) -> Counter:
    lines = text.splitlines()
    counts = Counter()
    for index, line in enumerate(lines):
        cells = _cells(line)
        if not cells or not all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells):
            continue
        width = len(cells)
        _require(index > 0 and len(_cells(lines[index - 1])) == width, f"{name}:{index + 1}: table separator/header column mismatch")
        counts["tables"] += 1
        counts["table_rows"] += 1  # Header, excluding separator.
        for row_index in range(index + 1, len(lines)):
            row = _cells(lines[row_index])
            if not row:
                break
            _require(len(row) == width, f"{name}:{row_index + 1}: table has {len(row)} columns; expected {width}")
            counts["table_rows"] += 1
    return counts


def _milestones(source: str, headings: list[tuple[int, int, str]]) -> dict:
    name = "learning/MILESTONES.md"
    records = [(line, title) for line, level, title in headings if level == 3 and re.match(r"M\d+\.", title)]
    numbers = [int(re.match(r"M(\d+)\.", title)[1]) for _, title in records]
    _require(numbers == list(range(1, 31)), f"{name}: expected M01 through M30 exactly once in order; got {numbers}")
    _require([title.split(".", 1)[0] for _, title in records] == [f"M{number:02}" for number in range(1, 31)], f"{name}: milestone identifiers must use M01 through M30")
    lines = source.splitlines()
    for index, (line_number, title) in enumerate(records):
        stop = records[index + 1][0] - 1 if index + 1 < len(records) else len(lines)
        block = "\n".join(lines[line_number:stop])
        for field in FIELDS:
            _require(block.count(f"**{field}**") == 1, f"{name}:{line_number}: {title} must have one {field} field")
        # The schema combines earlier dependencies and later story roles in one
        # paragraph. Only its first sentence declares the dependency list.
        declaration = block.split("**Dependencies and story role.**", 1)[1].split(".", 1)[0]
        dependencies = [int(value) for value in re.findall(r"\bM(\d{2})\b", declaration)]
        _require(all(1 <= dependency < index + 1 for dependency in dependencies), f"{name}:{line_number}: dependencies must precede M{index + 1:02}")
    split = {}
    for line_number, level, title in headings:
        if level != 2:
            continue
        interval = re.search(r"M(\d{2}) through M(\d{2})", title)
        if interval:
            first, last = map(int, interval.groups())
            next_sections = [line for line, other_level, _ in headings if other_level <= 2 and line > line_number]
            stop = min(next_sections, default=len(lines) + 1)
            members = [int(re.match(r"M(\d+)\.", value)[1]) for line, value in records if line_number < line < stop]
            _require(members == list(range(first, last + 1)), f"{name}:{line_number}: milestone section does not match its stated interval")
            split[f"M{first:02}-M{last:02}"] = len(members)
    _require(split == {"M01-M14": 14, "M15-M30": 16}, f"{name}: expected the structural 14/16 milestone split; got {split}")
    return {"count": len(records), "sections": split, "fields_per_milestone": len(FIELDS)}


def run_checks(root: Path | str | None = None) -> dict:
    root = Path(root).resolve() if root is not None else ROOT
    paths = sorted(path for path in root.rglob("*.md") if ".git" not in path.relative_to(root).parts)
    guide = root / "llms.txt"
    if guide.is_file():
        paths.append(guide)
    _require(bool(paths), f"{root}: no repository documents found")
    documents, headings, totals = {}, {}, Counter()
    for path in paths:
        name = path.relative_to(root).as_posix()
        source = path.read_text(encoding="utf-8")
        prose, math_masked, counts = _source_views(source, name)
        anchors, records = _headings(prose)
        documents[path] = (source, prose, anchors)
        headings[path] = records
        totals.update(counts)
        totals.update(_tables(math_masked, name))
    for path, (_, prose, _) in documents.items():
        name = path.relative_to(root).as_posix()
        for target, offset in _destinations(prose):
            url = urlsplit(target)
            if url.scheme or url.netloc:
                continue
            line_number = prose.count("\n", 0, offset) + 1
            destination = (path.parent / unquote(url.path)).resolve() if url.path else path
            _require(destination.is_relative_to(root), f"{name}:{line_number}: relative link leaves repository: {target}")
            _require(destination.exists(), f"{name}:{line_number}: missing link target: {target}")
            totals["relative_links"] += 1
            if url.fragment:
                _require(destination in documents, f"{name}:{line_number}: cannot verify heading fragment in {target}")
                anchor = unquote(url.fragment)
                _require(anchor in documents[destination][2], f"{name}:{line_number}: missing heading fragment: {target}")
                totals["heading_fragments"] += 1
    milestone_path = root / "learning/MILESTONES.md"
    _require(milestone_path in documents, "learning/MILESTONES.md: missing milestone map")
    milestone_report = _milestones(documents[milestone_path][0], headings[milestone_path])
    return {
        "status": "passed",
        "documents": len(documents),
        **dict(sorted(totals.items())),
        "milestones": milestone_report,
        "scope": "Source navigation and Markdown structure; hosted rendering requires visual inspection.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="Repository root to inspect.")
    parser.add_argument("--output", type=Path, help="Optional JSON report path.")
    args = parser.parse_args()
    report = run_checks(args.root)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
