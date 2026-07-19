#!/usr/bin/env python3
"""Deterministic IllinoisGRMHD knowledge-base checker."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import unicodedata
from collections import Counter, defaultdict, deque
from datetime import date as calendar_date
from pathlib import Path
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]
AGENTS = ROOT / "AGENTS.md"
WIKI = ROOT / "wiki"
SOURCES = ROOT / "raw/SOURCES.md"
CATALOG = WIKI / "catalog.md"
GLOSSARY = WIKI / "glossary.md"
SOURCE_MAP = WIKI / "source-map.md"
CONTRADICTIONS = WIKI / "contradictions.md"
ISOLATED_NAMESPACES = (WIKI / "grhaylhd", WIKI / "grhaylid")

LEAF_SECTIONS = ["Summary", "Detail", "Sources", "See Also"]
PAGE_STATUSES = {"confirmed", "provisional", "contested", "stale"}
AUTHORITY_TIERS = {
    "primary-code",
    "primary-test",
    "generated-evidence",
    "primary-doc",
    "background",
}
COVERAGE_STATES = {"unmapped", "partial", "complete"}
GLOBAL_WIKI = {
    "SCHEMA.md",
    "catalog.md",
    "contradictions.md",
    "glossary.md",
    "source-map.md",
    "workflows.md",
    "lint/CHECKS.md",
}
SUPPORT_PAGES = {WIKI / path for path in GLOBAL_WIKI}
BRANCH_CHILDREN = {
    "architecture": {
        "overview.md",
        "cactus-surface-and-build.md",
        "schedule-lifecycle.md",
    },
    "evolution": {
        "state-and-eos-modes.md",
        "primitive-conservative-conversion.md",
        "con2prim-recovery-and-diagnostics.md",
        "reconstruction-fluxes-and-sources.md",
        "matter-boundaries-and-perturbations.md",
    },
    "magnetics": {
        "staggered-state-and-magnetic-reconstruction.md",
        "induction-and-lorenz-gauge-rhs.md",
        "electromagnetic-boundaries-and-symmetry.md",
    },
    "integration": {
        "hydrobase-grhaylib-and-tmunu.md",
        "parameters-and-runtime-controls.md",
        "migration-and-backward-compatibility.md",
    },
    "validation": {
        "test-harness-and-oracles.md",
        "balsara-and-tov-cases.md",
    },
}
EXPECTED_WIKI = GLOBAL_WIKI | {
    path
    for branch, children in BRANCH_CHILDREN.items()
    for path in ({f"{branch}/index.md"} | {f"{branch}/{child}" for child in children})
}
ROUTER_CHILDREN = {
    "AGENTS.md": {f"wiki/{branch}/index.md" for branch in BRANCH_CHILDREN},
    **{
        f"wiki/{branch}/index.md": {
            f"wiki/{branch}/{child}" for child in children
        }
        for branch, children in BRANCH_CHILDREN.items()
    },
}
CATALOG_COLUMNS = [
    "Page",
    "Type",
    "One-line answer",
    "Route",
    "Tags",
    "Aliases / Query terms",
    "Status",
    "Last reconciled",
]
SOURCE_MAP_COLUMNS = [
    "Source / aggregate",
    "Authority tier",
    "Dependency coverage",
    "Dependent pages",
    "Covered paths / stable locators",
    "Known gaps",
    "Last check",
    "Next action",
]
CONTRADICTION_COLUMNS = [
    "ID",
    "Claim",
    "Claim status",
    "Source A",
    "Source B",
    "Authority decision",
    "Affected pages",
    "Page-status rationale",
    "Owner/trigger",
    "Resolution test",
    "Opened",
    "Resolved",
    "Notes",
]
COMMISSIONED = {
    "plan1.md",
    "plan2.md",
    "plan3.md",
    "plan_synth.md",
    "tasks1.md",
    "tasks2.md",
    "tasks3.md",
    "tasks4.md",
}
LINK_RE = re.compile(r"(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)")
SOURCE_PATH_RE = re.compile(
    r"(?<![\w/.])(?:(?:\.\./)+|\./)?"
    r"(?P<source>IllinoisGRMHD/[A-Za-z0-9_./*?+\[\]-]+)"
)
MANIFEST_ID_RE = re.compile(
    r"(?<![\w-])(illinoisgrmhd-[a-z0-9]+(?:-[a-z0-9]+)*)(?![\w-])"
)
INLINE_LITERAL_RE = re.compile(r"`([^`\n]+)`")
RETAINED_DATE_RE = re.compile(
    r"\bLast (?:reconciled|checked|audited):[ \t]*([^ \t·.]*)"
)
FINGERPRINT_VALUE_RE = re.compile(
    r"\b(?:sha-?\d+|md5|blake\w*|hash|digest)\b"
    r"(?:\s*[:=]\s*|\s+)[`*_~]{0,3}(?:0x)?[0-9a-f]{8,}\b",
    re.IGNORECASE,
)


def rel(path: Path) -> str:
    """Return repository-relative POSIX path when possible."""
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def inside(path: Path, parent: Path) -> bool:
    """Return whether resolved path lies at or below parent."""
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError:
        return False
    return True


def in_isolated(path: Path) -> bool:
    """Return whether path belongs to an additive evidence namespace."""
    return any(inside(path, namespace) for namespace in ISOLATED_NAMESPACES)


def fail(failures: list[str], path: Path, message: str, line: int | None = None) -> None:
    """Append stable diagnostic."""
    location = rel(path)
    if line is not None:
        location += f":{line}"
    failures.append(f"{location}: {message}")


def read(path: Path) -> str:
    """Read UTF-8 text."""
    return path.read_text(encoding="utf-8")


def governed_files() -> list[Path]:
    """Return governed text files."""
    files = [AGENTS, SOURCES, Path(__file__).resolve()]
    grhaylhd_sources = ROOT / "raw/grhaylhd/SOURCES.md"
    if grhaylhd_sources.is_file():
        files.append(grhaylhd_sources)
    grhaylid_sources = ROOT / "raw/grhaylid/SOURCES.md"
    if grhaylid_sources.is_file():
        files.append(grhaylid_sources)
    files.extend(sorted(WIKI.rglob("*.md")))
    return sorted({path.resolve() for path in files})


def wiki_pages() -> list[Path]:
    """Return live wiki Markdown pages."""
    return sorted(path.resolve() for path in WIKI.rglob("*.md"))


def mask_code(text: str) -> str:
    """Mask fenced code and HTML comments while preserving line count."""
    def masked(value: str) -> str:
        return "".join("\n" if char == "\n" else " " for char in value)

    html_comment = re.compile(r"<!--.*?(?:-->|$)", re.DOTALL)
    text = html_comment.sub(lambda match: masked(match.group(0)), text)

    result: list[str] = []
    fence_char: str | None = None
    fence_length = 0
    for line in text.splitlines(keepends=True):
        content = line.rstrip("\r\n")
        if fence_char is not None:
            result.append(masked(line))
            closer = re.fullmatch(
                rf"[ ]{{0,3}}{re.escape(fence_char)}"
                rf"{{{fence_length},}}[ \t]*",
                content,
            )
            if closer is not None:
                fence_char = None
                fence_length = 0
            continue

        opener = re.match(r"^[ ]{0,3}(`{3,}|~{3,})(.*)$", content)
        if opener is not None:
            marker, info = opener.groups()
            if marker[0] != "`" or "`" not in info:
                fence_char = marker[0]
                fence_length = len(marker)
                result.append(masked(line))
                continue
        result.append(line)
    return "".join(result)


def mask_inline_code(text: str) -> str:
    """Mask inline-code spans while preserving offsets."""
    pattern = re.compile(
        r"(?<!`)(?P<ticks>`+)(?!`)(?P<body>.*?)(?<!`)(?P=ticks)(?!`)",
        re.DOTALL,
    )

    def repl(match: re.Match[str]) -> str:
        return "".join("\n" if char == "\n" else " " for char in match.group(0))

    return pattern.sub(repl, text)


def slug(title: str) -> str:
    """Return GitHub-style heading slug sufficient for this KB."""
    value = re.sub(r"[^\w\s-]", "", title.strip().lower())
    return re.sub(r"\s+", "-", value)


def anchors(path: Path) -> set[str]:
    """Return Markdown heading anchors."""
    values: set[str] = set()
    counts: Counter[str] = Counter()
    for line in mask_code(read(path)).splitlines():
        match = re.match(r"^#{1,6}\s+(.+?)\s*$", line)
        if not match:
            continue
        base = slug(re.sub(r"\s+#+\s*$", "", match.group(1)))
        number = counts[base]
        counts[base] += 1
        values.add(base if number == 0 else f"{base}-{number}")
    return values


def resolve_link(source: Path, target: str) -> tuple[Path | None, str]:
    """Resolve local Markdown link target and anchor."""
    parsed = urlparse(target.strip())
    if parsed.scheme in {"http", "https", "mailto", "tel"} or parsed.netloc:
        return None, ""
    raw_path = unquote(parsed.path)
    anchor = unquote(parsed.fragment)
    if not raw_path:
        return source.resolve(), anchor
    if raw_path.startswith("/"):
        return (ROOT / raw_path.lstrip("/")).resolve(), anchor
    return (source.parent / raw_path).resolve(), anchor


def links(path: Path) -> list[tuple[int, str, str, Path | None, str]]:
    """Return ordinary single-line Markdown links."""
    text = mask_inline_code(mask_code(read(path)))
    result = []
    for match in LINK_RE.finditer(text):
        target, anchor = resolve_link(path, match.group(2))
        result.append(
            (
                text.count("\n", 0, match.start()) + 1,
                match.group(1),
                match.group(2),
                target,
                anchor,
            )
        )
    return result


def exact_link_match(cell: str) -> re.Match[str] | None:
    """Return one ordinary Markdown link only when it fills the table cell."""
    return LINK_RE.fullmatch(mask_inline_code(cell))


def source_path_mentions(text: str) -> list[tuple[int, str]]:
    """Return IllinoisGRMHD path/glob mentions independent of Markdown style."""
    masked = mask_code(text)
    return [
        (
            masked.count("\n", 0, match.start()) + 1,
            match.group("source").rstrip(".,;:"),
        )
        for match in SOURCE_PATH_RE.finditer(masked)
    ]


def manifest_id_mentions(text: str) -> list[tuple[int, str]]:
    """Return manifest-ID mentions independent of Markdown style."""
    masked = mask_code(text)
    return [
        (masked.count("\n", 0, match.start()) + 1, match.group(1))
        for match in MANIFEST_ID_RE.finditer(masked)
    ]


def inline_literals(text: str) -> list[tuple[int, str]]:
    """Return inline-code literals outside fenced code and HTML comments."""
    masked = mask_code(text)
    return [
        (masked.count("\n", 0, match.start()) + 1, match.group(1))
        for match in INLINE_LITERAL_RE.finditer(masked)
    ]


def tables(path: Path) -> list[tuple[list[str], list[tuple[int, list[str]]]]]:
    """Parse simple Markdown tables with a required delimiter row."""
    result: list[tuple[list[str], list[tuple[int, list[str]]]]] = []
    header: list[str] | None = None
    rows: list[tuple[int, list[str]]] = []
    separator_seen = False

    for number, line in enumerate(mask_code(read(path)).splitlines(), start=1):
        if not line.startswith("|"):
            if header is not None and separator_seen:
                result.append((header, rows))
            header, rows, separator_seen = None, [], False
            continue

        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        is_separator = bool(cells) and all(
            re.fullmatch(r":?-{3,}:?", cell) for cell in cells
        )
        if header is None:
            header = cells
        elif not separator_seen:
            if is_separator and len(cells) == len(header):
                separator_seen = True
            else:
                header, rows, separator_seen = None, [], False
        elif not is_separator:
            rows.append((number, cells))

    if header is not None and separator_seen:
        result.append((header, rows))
    return result


def table(path: Path, first_column: str) -> tuple[list[str], list[tuple[int, list[str]]]]:
    """Return table with requested first-column label."""
    for header, rows in tables(path):
        if header and header[0] == first_column:
            return header, rows
    return [], []


def named_section_table(
    path: Path, heading: str, expected_header: list[str]
) -> tuple[bool, list[tuple[int, list[str]]]]:
    """Return one exact simple table from one uniquely named H2 section."""
    lines = mask_code(read(path)).splitlines()
    headings = [
        index
        for index, line in enumerate(lines)
        if line == f"## {heading}"
    ]
    if len(headings) != 1:
        return False, []
    start = headings[0] + 1
    end = next(
        (
            index
            for index in range(start, len(lines))
            if re.match(r"^##\s+", lines[index])
        ),
        len(lines),
    )

    def cells(line: str) -> list[str]:
        if not line.startswith("|"):
            return []
        return [cell.strip() for cell in line.strip().strip("|").split("|")]

    headers = [
        index
        for index in range(start, end)
        if cells(lines[index]) == expected_header
    ]
    if len(headers) != 1 or headers[0] + 1 >= end:
        return False, []
    separator = cells(lines[headers[0] + 1])
    if len(separator) != len(expected_header) or not all(
        re.fullmatch(r":?-{3,}:?", cell) for cell in separator
    ):
        return False, []
    rows: list[tuple[int, list[str]]] = []
    for index in range(headers[0] + 2, end):
        row = cells(lines[index])
        if not row:
            break
        rows.append((index + 1, row))
    return True, rows


def valid_date(value: str, allow_na: bool = False) -> bool:
    """Validate MM-DD-YYYY calendar date."""
    if allow_na and value.lower() in {"n/a", "-"}:
        return True
    match = re.fullmatch(r"(\d{2})-(\d{2})-(\d{4})", value)
    if not match:
        return False
    month, day, year = map(int, match.groups())
    try:
        calendar_date(year, month, day)
    except ValueError:
        return False
    return True


def page_metadata(path: Path) -> tuple[str | None, str | None]:
    """Extract page status/date from first five lines."""
    header = "\n".join(read(path).splitlines()[:5])
    status = re.search(r"\bStatus:\s*(\w+)\b", header)
    dated = re.search(r"\bLast (?:reconciled|checked|audited):\s*([^\s·]+)", header)
    return (status.group(1).lower() if status else None, dated.group(1) if dated else None)


def section(text: str, name: str) -> str:
    """Return H2 section body."""
    match = re.search(rf"^## {re.escape(name)}\s*$", text, re.MULTILINE)
    if not match:
        return ""
    next_heading = re.search(r"^##\s+", text[match.end() :], re.MULTILINE)
    end = match.end() + next_heading.start() if next_heading else len(text)
    return text[match.end() : end].strip()


def check_inventory(failures: list[str]) -> None:
    """Check exact commissioned KB inventory."""
    for path in (AGENTS, SOURCES, Path(__file__).resolve()):
        if not path.is_file():
            fail(failures, path, "required KB node missing")
    live = {
        path.relative_to(WIKI).as_posix()
        for path in wiki_pages()
        if not in_isolated(path)
    }
    for missing in sorted(EXPECTED_WIKI - live):
        fail(failures, WIKI / missing, "required wiki page missing")
    for extra in sorted(live - EXPECTED_WIKI):
        fail(failures, WIKI / extra, "unexpected wiki page outside exact tree")
    raw_files = {
        path.relative_to(ROOT / "raw").as_posix()
        for path in (ROOT / "raw").rglob("*")
        if path.is_file()
    }
    for extra in sorted(
        path for path in raw_files - {"SOURCES.md"}
        if not path.startswith(("grhaylhd/", "grhaylid/"))
    ):
        fail(failures, ROOT / "raw" / extra, "unexpected file outside exact KB tree")


def check_text(failures: list[str]) -> None:
    """Check governed text hygiene and prohibited artifacts."""
    for path in governed_files():
        data = path.read_bytes()
        if not data or not data.endswith(b"\n"):
            fail(failures, path, "missing final newline")
        text = data.decode("utf-8")
        masked_lines = mask_code(text).splitlines()
        for number, line in enumerate(text.splitlines(), start=1):
            if re.search(r"[ \t]+$", line):
                fail(failures, path, "trailing horizontal whitespace", number)
            if path.suffix != ".md":
                continue
            masked_line = masked_lines[number - 1]
            if re.search(r"\[\[[^\]]+\]\]", masked_line):
                fail(failures, path, "Obsidian wikilink found", number)
            if re.search(r"\b(?:TODO|TBD|PLACEHOLDER|XX-XX-XXXX)\b", line, re.IGNORECASE):
                fail(failures, path, "unfilled template marker", number)
            if re.search(r"\b\d{4}-\d{2}-\d{2}\b", line):
                fail(failures, path, "YYYY-MM-DD date literal found", number)
            for retained in RETAINED_DATE_RE.finditer(line):
                if not valid_date(retained.group(1)):
                    fail(
                        failures,
                        path,
                        f"invalid retained date: {retained.group(1)!r}",
                        number,
                    )
            for candidate in re.findall(r"\b\d{2}-\d{2}-\d{4}\b", line):
                if not valid_date(candidate):
                    fail(failures, path, f"invalid MM-DD-YYYY date: {candidate!r}", number)
            if FINGERPRINT_VALUE_RE.search(line):
                fail(failures, path, "stored source fingerprint found", number)
            if re.search(r"\|\s*(?:Mtime|Hash|Digest)\s*\|", line, re.IGNORECASE):
                fail(failures, path, "source fingerprint/modification-time column found", number)
            if re.search(r"\b(?:mtime|modification[- ]time)\s*[:=]\s*\S+", line, re.IGNORECASE):
                fail(failures, path, "stored source modification-time value found", number)
    for base in (WIKI, ROOT / "raw"):
        for path in sorted(base.rglob("*")):
            if path.is_file() and path.suffix.lower() != ".md":
                fail(failures, path, "non-Markdown artifact in KB tree")
    for path in sorted(ROOT.glob("*.md")):
        if path.name in {"AGENTS.md", "README.md", "CITATION.md"} | COMMISSIONED:
            continue
        stem = path.stem.casefold()
        is_plan_or_task = bool(
            re.fullmatch(r"(?:plan|tasks?)(?:\d+|[._-].*)?", stem)
        )
        is_other_maintenance = bool(
            re.search(
                r"(?:^|[_-])(?:audit|report|scan|rank|log)(?:[_.-]|$)",
                path.name,
                re.IGNORECASE,
            )
        )
        if is_plan_or_task or is_other_maintenance:
            fail(failures, path, "uncommissioned root maintenance artifact")


def check_links(failures: list[str]) -> dict[Path, set[Path]]:
    """Check local links and return wiki reachability graph."""
    graph: dict[Path, set[Path]] = defaultdict(set)
    anchor_cache: dict[Path, set[str]] = {}
    for path in governed_files():
        if path.suffix != ".md":
            continue
        for number, _, raw, target, anchor in links(path):
            if target is None:
                continue
            if not inside(target, ROOT):
                fail(failures, path, f"link target escapes repository: {raw}", number)
                continue
            if not target.exists():
                fail(failures, path, f"link target missing: {raw}", number)
                continue
            if anchor:
                if not target.is_file():
                    fail(failures, path, f"link anchor target is not a file: {raw}", number)
                else:
                    anchor_cache.setdefault(target, anchors(target))
                    if slug(anchor) not in anchor_cache[target]:
                        fail(failures, path, f"link anchor missing: {raw}", number)
                        continue
            if target.is_file() and (target == AGENTS or inside(target, WIKI)):
                graph[path.resolve()].add(target.resolve())
            if target.name in COMMISSIONED:
                fail(failures, path, f"coordination artifact enters KB routing: {target.name}", number)
    return graph


def check_routers(failures: list[str]) -> None:
    """Check exact router shape and immediate-child edges."""
    allowed_global = {path.resolve() for path in SUPPORT_PAGES} | {
        AGENTS.resolve(), SOURCES.resolve()
    }
    for router_name, child_names in ROUTER_CHILDREN.items():
        router = (ROOT / router_name).resolve()
        if not router.exists():
            continue
        text = mask_code(read(router))
        h1 = re.findall(r"^#\s+.+$", text, re.MULTILINE)
        if len(h1) != 1:
            fail(failures, router, "router must contain exactly one H1")
        if router != AGENTS:
            quote_lines = [line for line in text.splitlines() if line.startswith(">")]
            if len(quote_lines) != 1:
                fail(failures, router, "router must contain one-line blockquote")
            header, _ = table(router, "Page")
            if header != ["Page", "Go here when..."]:
                fail(failures, router, "router lacks exact Page / Go here when... table")
        if re.search(r"^## (?:Detail|Sources)\s*$", text, re.MULTILINE):
            fail(failures, router, "router has Detail or Sources section")
        child_paths = {(ROOT / child).resolve() for child in child_names}
        targets = {
            target
            for _, _, _, target, _ in links(router)
            if target is not None and inside(target, ROOT)
        }
        allowed_for_router = allowed_global | (
            {
                (GR_WIKI / item).resolve()
                for item in GR_ROOT_ROUTES
            } | {
                (GID_WIKI / item).resolve()
                for item in GID_ROOT_ROUTES
            } | {GR_REGISTRY.resolve(), GID_REGISTRY.resolve()}
            if router == AGENTS.resolve() else set()
        )
        for missing in sorted(child_paths - targets):
            fail(failures, router, f"router missing immediate child: {rel(missing)}")
        for target in sorted(targets):
            if target not in child_paths and target not in allowed_for_router:
                fail(failures, router, f"illegal router edge: {rel(target)}")


def check_leaves(failures: list[str]) -> None:
    """Check leaf header, section, parent, source, and neighbor contract."""
    for path in wiki_pages():
        if in_isolated(path):
            continue
        if path.name == "index.md" or path in {item.resolve() for item in SUPPORT_PAGES}:
            continue
        text = mask_code(read(path))
        lines = text.splitlines()
        h1_lines = [line for line in lines if re.fullmatch(r"#\s+.+", line)]
        status_match = re.fullmatch(
            r"> Status: (confirmed|provisional|contested|stale) · "
            r"Last reconciled: (\d{2}-\d{2}-\d{4})",
            lines[2] if len(lines) > 2 else "",
        )
        up_match = re.fullmatch(
            r"> Up: \[[^\]\n]+\]\([^)\n]+\)",
            lines[3] if len(lines) > 3 else "",
        )
        boundary_ok = (
            len(lines) > 5
            and len(h1_lines) == 1
            and lines[0] == h1_lines[0]
            and lines[1] == ""
            and lines[4] == ""
            and lines[5] == "## Summary"
        )
        if not boundary_ok:
            fail(
                failures,
                path,
                "leaf must begin H1, blank, Status, Up, blank, then immediate Summary",
            )
        if status_match is None or not valid_date(status_match.group(2)):
            fail(failures, path, "leaf line 3 must be exact Status/date header")
        if up_match is None:
            fail(failures, path, "leaf line 4 must be exact Up header")
        quote_lines = lines[2:4]
        headings = [match.group(1).strip() for match in re.finditer(r"^##\s+(.+?)\s*$", text, re.MULTILINE)]
        if headings != LEAF_SECTIONS:
            fail(failures, path, "leaf H2 sections must be exactly Summary, Detail, Sources, See Also")
        for name in LEAF_SECTIONS:
            if not section(text, name):
                fail(failures, path, f"leaf has empty {name}")
        parent = (path.parent / "index.md").resolve()
        up_targets: list[Path] = []
        for line in quote_lines:
            if line.startswith("> Up:"):
                for match in LINK_RE.finditer(line):
                    target, _ = resolve_link(path, match.group(2))
                    if target is not None:
                        up_targets.append(target)
        if up_targets != [parent]:
            fail(failures, path, f"leaf Up link must resolve to {rel(parent)}")
        see_also = section(text, "See Also")
        see_targets = {
            target
            for match in LINK_RE.finditer(see_also)
            for target, _ in [resolve_link(path, match.group(2))]
            if target is not None and inside(target, WIKI)
        }
        if parent not in see_targets:
            fail(failures, path, "leaf See Also lacks parent router")
        live = set(wiki_pages())
        typed_targets: set[Path] = set()
        typed_prefix = re.compile(
            r"^-\s+(?:Depends on|Implements|Example|Contrasts with|See also):\s+",
            re.IGNORECASE,
        )
        for line in see_also.splitlines():
            prefix = typed_prefix.match(line)
            if prefix is None:
                continue
            for match in LINK_RE.finditer(line[prefix.end() :]):
                target, _ = resolve_link(path, match.group(2))
                if target in live and target not in {path.resolve(), parent}:
                    typed_targets.add(target)
        if not typed_targets:
            fail(failures, path, "leaf See Also lacks typed neighboring wiki page")


def parse_manifest(failures: list[str]) -> tuple[dict[str, str], set[str]]:
    """Parse manifest, validate patterns, and return ID/pattern plus covered files."""
    header, rows = table(SOURCES, "Source")
    expected = ["Source", "Provenance", "Status", "Ingest"]
    if header != expected:
        fail(failures, SOURCES, "manifest columns do not match schema")
        return {}, set()
    identifiers: dict[str, str] = {}
    provenances: set[str] = set()
    covered: set[str] = set()
    for number, row in rows:
        if len(row) != 4:
            fail(failures, SOURCES, "manifest row has wrong field count", number)
            continue
        if any(not cell for cell in row):
            fail(failures, SOURCES, "manifest row contains empty cell", number)
            continue
        source_match = re.fullmatch(r"`([^`\n]+)`", row[0])
        provenance_match = re.fullmatch(r"`([^`\n]+)`", row[1])
        if not source_match or not provenance_match or row[1].count("`") != 2:
            fail(failures, SOURCES, "manifest Source/Provenance must each contain one literal", number)
            continue
        identifier, pattern = source_match.group(1), provenance_match.group(1)
        if identifier in identifiers:
            fail(failures, SOURCES, f"duplicate manifest source ID: {identifier}", number)
        if pattern in provenances:
            fail(failures, SOURCES, f"duplicate manifest provenance: {pattern}", number)
        identifiers[identifier] = pattern
        provenances.add(pattern)
        if row[2] not in {"living", "frozen"}:
            fail(failures, SOURCES, f"invalid manifest Status: {row[2]!r}", number)
        if row[3] not in {"registered", "partial", "ingested"}:
            fail(failures, SOURCES, f"invalid manifest Ingest: {row[3]!r}", number)
        components = pattern.split("/")
        invalid = (
            not pattern
            or "\\" in pattern
            or pattern.startswith("/")
            or bool(urlparse(pattern).scheme)
            or any(component in {"", ".", ".."} for component in components)
            or components[0] != "IllinoisGRMHD"
            or any("**" in component and component != "**" for component in components)
        )
        if invalid:
            fail(failures, SOURCES, f"invalid manifest provenance: {pattern}", number)
            continue
        try:
            matches = sorted(ROOT.glob(pattern), key=lambda item: item.as_posix())
        except (OSError, RuntimeError, ValueError, NotImplementedError) as error:
            fail(failures, SOURCES, f"invalid manifest provenance: {pattern}: {error}", number)
            continue
        if not matches:
            fail(failures, SOURCES, f"manifest pattern unmatched: {pattern}", number)
            continue
        files = [item for item in matches if item.is_file()]
        if not files:
            fail(failures, SOURCES, f"manifest pattern matches directories only: {pattern}", number)
            continue
        for item in files:
            if not inside(item, ROOT) or not inside(item, ROOT / "IllinoisGRMHD"):
                fail(failures, SOURCES, f"manifest match escapes IllinoisGRMHD: {rel(item)}", number)
                continue
            covered.add(rel(item))
    command = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "--", "IllinoisGRMHD"],
        text=True,
        capture_output=True,
        check=False,
    )
    if command.returncode != 0:
        fail(failures, ROOT / ".git", f"git ls-files failed: {command.stderr.strip()}")
        tracked_files: set[str] = set()
    else:
        tracked_files = set(filter(None, command.stdout.splitlines()))
        for tracked in sorted(tracked_files):
            if tracked not in covered:
                fail(failures, ROOT / tracked, "tracked IllinoisGRMHD file is uncovered")
    return identifiers, covered & tracked_files


def source_section_for(path: Path) -> str:
    """Return original Sources body."""
    return section(mask_code(read(path)), "Sources")


def literal_local_target(source: Path, literal: str) -> Path | None:
    """Resolve an existing path-like inline literal using KB path conventions."""
    candidate = literal.split("::", 1)[0].strip().rstrip(".,;:")
    if "/" not in candidate:
        return None
    if candidate.startswith("IllinoisGRMHD/"):
        target = ROOT / candidate
    elif candidate.startswith("/"):
        target = Path(candidate)
    elif candidate.startswith("."):
        target = source.parent / candidate
    else:
        target = ROOT / candidate
    try:
        resolved = target.resolve()
        return resolved if resolved.exists() else None
    except (OSError, RuntimeError, ValueError):
        return None


def registered_source_mention(candidate: str, covered: set[str]) -> bool:
    """Return whether a path/glob mention reaches tracked manifest coverage."""
    if any(char in candidate for char in "*?["):
        try:
            return any(
                rel(item) in covered
                for item in ROOT.glob(candidate)
                if item.is_file()
            )
        except (OSError, RuntimeError, ValueError, NotImplementedError):
            return False
    return candidate in covered or any(
        item.startswith(candidate.rstrip("/") + "/") for item in covered
    )


def check_source_citations(
    failures: list[str], identifiers: dict[str, str], covered: set[str]
) -> None:
    """Check all leaf local citations stay in scope and Sources has evidence."""
    allowed_support = {AGENTS.resolve(), SOURCES.resolve(), Path(__file__).resolve()}
    for path in wiki_pages():
        if in_isolated(path):
            continue
        if path.name == "index.md" or path in {item.resolve() for item in SUPPORT_PAGES}:
            continue

        text = read(path)
        for number, identifier in manifest_id_mentions(text):
            if identifier not in identifiers:
                fail(
                    failures,
                    path,
                    f"unknown manifest source ID: {identifier}",
                    number,
                )

        for number, _, raw, target, _ in links(path):
            if target is None or inside(target, WIKI):
                continue
            if target in allowed_support:
                continue
            if not inside(target, ROOT / "IllinoisGRMHD"):
                fail(failures, path, f"local citation outside IllinoisGRMHD: {raw}", number)

        for number, literal in inline_literals(text):
            target = literal_local_target(path, literal)
            if target is None or inside(target, WIKI) or target in allowed_support:
                continue
            if not inside(target, ROOT / "IllinoisGRMHD"):
                fail(
                    failures,
                    path,
                    f"local literal citation outside IllinoisGRMHD: {literal}",
                    number,
                )

        body = source_section_for(path)
        evidence = False
        for match in LINK_RE.finditer(body):
            parsed = urlparse(match.group(2))
            if parsed.scheme in {"http", "https"}:
                continue
            target, _ = resolve_link(path, match.group(2))
            if target is None:
                continue
            if not inside(target, ROOT / "IllinoisGRMHD"):
                fail(failures, path, f"source outside IllinoisGRMHD: {match.group(2)}")
                continue
            literal = rel(target)
            registered = literal in covered or (
                target.is_dir() and any(item.startswith(literal.rstrip("/") + "/") for item in covered)
            )
            if not registered:
                fail(failures, path, f"source is not registered: {literal}")
            else:
                evidence = True
        for _, candidate in source_path_mentions(body):
            if registered_source_mention(candidate, covered):
                evidence = True
            else:
                fail(failures, path, f"source is not registered: {candidate}")
        for _, identifier in manifest_id_mentions(body):
            if identifier in identifiers:
                evidence = True
        if not evidence:
            fail(failures, path, "leaf Sources has no registered IllinoisGRMHD evidence")


def check_reachability(failures: list[str], graph: dict[Path, set[Path]]) -> None:
    """Check every wiki page is reachable from root."""
    reached: set[Path] = set()
    queue = deque([AGENTS.resolve()])
    while queue:
        current = queue.popleft()
        if current in reached:
            continue
        reached.add(current)
        queue.extend(sorted(graph.get(current, set()) - reached))
    for path in wiki_pages():
        if in_isolated(path):
            continue
        if path not in reached:
            fail(failures, path, "wiki page not reachable from AGENTS.md")


def check_catalog(failures: list[str]) -> None:
    """Check catalog exact live set and page metadata agreement."""
    header, rows = table(CATALOG, "Page")
    if header != CATALOG_COLUMNS:
        fail(failures, CATALOG, "catalog columns do not match schema")
        return
    seen: list[tuple[Path, int, list[str]]] = []
    for number, row in rows:
        if len(row) != len(header):
            fail(failures, CATALOG, "catalog row has wrong field count", number)
            continue
        matches = list(LINK_RE.finditer(row[0]))
        if len(matches) != 1:
            fail(failures, CATALOG, "catalog Page cell must contain one link", number)
            continue
        target, _ = resolve_link(CATALOG, matches[0].group(2))
        if target is None or not inside(target, WIKI):
            fail(failures, CATALOG, "catalog Page target is outside wiki", number)
            continue
        seen.append((target, number, row))
    counts = Counter(target for target, _, _ in seen)
    live = {path for path in wiki_pages() if not in_isolated(path)}
    for path in sorted(live):
        if counts[path] == 0:
            fail(failures, CATALOG, f"missing wiki page: {rel(path)}")
        elif counts[path] > 1:
            fail(failures, CATALOG, f"duplicate wiki page: {rel(path)}")
    for target, number, row in seen:
        if target not in live:
            fail(failures, CATALOG, f"catalog target is not a live wiki page: {rel(target)}", number)
            continue
        if target.name == "index.md":
            expected_type, expected_status, expected_date = "router", "router", "n/a"
        else:
            expected_type = "leaf" if target not in {item.resolve() for item in SUPPORT_PAGES} else "governance"
            expected_status, expected_date = page_metadata(target)
        if row[1].lower() != expected_type:
            fail(failures, CATALOG, f"catalog type disagrees with {rel(target)}", number)
        if row[6].lower() != expected_status:
            fail(failures, CATALOG, f"catalog status disagrees with {rel(target)}", number)
        if row[7].lower() != str(expected_date).lower():
            fail(failures, CATALOG, f"catalog date disagrees with {rel(target)}", number)


def direct_source_dependents(
    identifiers: dict[str, str], covered: set[str]
) -> dict[str, set[Path]]:
    """Map IDs to pages directly citing validated, tracked manifest files."""
    glob_errors = (OSError, RuntimeError, ValueError, NotImplementedError)
    covered_paths = {(ROOT / item).resolve() for item in covered}

    def validated_expansion(pattern: str) -> set[Path]:
        try:
            return {
                item.resolve()
                for item in ROOT.glob(pattern)
                if item.is_file() and item.resolve() in covered_paths
            }
        except glob_errors:
            return set()

    def covered_targets(target: Path) -> set[Path]:
        try:
            resolved = target.resolve()
            if resolved in covered_paths:
                return {resolved}
            if target.is_dir():
                return {item for item in covered_paths if inside(item, resolved)}
        except glob_errors:
            pass
        return set()

    source_files = {
        identifier: validated_expansion(pattern)
        for identifier, pattern in identifiers.items()
    }
    expected = {identifier: set() for identifier in identifiers}
    support = {item.resolve() for item in SUPPORT_PAGES}

    for path in wiki_pages():
        if in_isolated(path):
            continue
        if path.name == "index.md" or (path in support and path != CONTRADICTIONS.resolve()):
            continue
        cited: set[Path] = set()
        for _, _, _, target, _ in links(path):
            if target is None or not inside(target, ROOT / "IllinoisGRMHD"):
                continue
            cited.update(covered_targets(target))

        text = read(path)
        for _, identifier in manifest_id_mentions(text):
            if identifier in source_files and source_files[identifier]:
                expected[identifier].add(path)

        for _, candidate in source_path_mentions(text):
            try:
                matches = (
                    ROOT.glob(candidate)
                    if any(char in candidate for char in "*?[")
                    else [ROOT / candidate]
                )
                for item in matches:
                    cited.update(covered_targets(item))
            except glob_errors:
                continue

        for identifier, files in source_files.items():
            if files & cited:
                expected[identifier].add(path)
    return expected


def check_source_map(
    failures: list[str], identifiers: dict[str, str], covered: set[str]
) -> None:
    """Check source-map identity, vocabulary, target, and coverage contract."""
    header, rows = table(SOURCE_MAP, "Source / aggregate")
    if header != SOURCE_MAP_COLUMNS:
        fail(failures, SOURCE_MAP, "source-map columns do not match schema")
        return
    prose = read(SOURCE_MAP)
    if re.search(r"\b(?:registered|ingested)\b", prose, re.IGNORECASE) or any(
        cell.lower() == "ingest" for cell in header
    ):
        fail(failures, SOURCE_MAP, "source map contains ingest-state vocabulary")
    expected_dependents = direct_source_dependents(identifiers, covered)
    seen: Counter[str] = Counter()
    live = set(wiki_pages())
    for number, row in rows:
        if len(row) != len(header):
            fail(failures, SOURCE_MAP, "source-map row has wrong field count", number)
            continue
        match = re.fullmatch(r"`([^`\n]+)`", row[0])
        if not match:
            fail(failures, SOURCE_MAP, "source-map source cell must contain one literal", number)
            continue
        identifier = match.group(1)
        seen[identifier] += 1
        if identifier not in identifiers:
            fail(failures, SOURCE_MAP, f"source map has unknown manifest ID: {identifier}", number)
        tiers = {value.strip() for value in row[1].split(";")}
        if not tiers or not tiers <= AUTHORITY_TIERS:
            fail(failures, SOURCE_MAP, "invalid source authority tier", number)
        coverage = row[2].lower()
        if coverage not in COVERAGE_STATES:
            fail(failures, SOURCE_MAP, f"invalid dependency coverage: {row[2]!r}", number)
        dependent_targets: set[Path] = set()
        for link in LINK_RE.finditer(row[3]):
            target, _ = resolve_link(SOURCE_MAP, link.group(2))
            if target is None or target not in live:
                fail(failures, SOURCE_MAP, "source-map dependent target is not live wiki page", number)
            else:
                dependent_targets.add(target)
        if coverage == "complete" and not dependent_targets:
            fail(failures, SOURCE_MAP, "complete source-map row lacks valid dependent pages", number)
        if coverage == "complete":
            for missing in sorted(expected_dependents.get(identifier, set()) - dependent_targets):
                fail(
                    failures,
                    SOURCE_MAP,
                    f"complete source-map row misses direct dependent: {rel(missing)}",
                    number,
                )
        if not valid_date(row[6]):
            fail(failures, SOURCE_MAP, f"invalid source-map Last check: {row[6]!r}", number)
    for identifier in sorted(identifiers):
        if seen[identifier] == 0:
            fail(failures, SOURCE_MAP, f"source map missing manifest ID: {identifier}")
        elif seen[identifier] > 1:
            fail(failures, SOURCE_MAP, f"source map duplicates manifest ID: {identifier}")


def check_glossary(failures: list[str]) -> None:
    """Check one unique live owner per glossary term."""
    header, rows = table(GLOSSARY, "Term")
    if header != ["Term", "Routing meaning", "Owner"]:
        fail(failures, GLOSSARY, "glossary columns do not match schema")
        return
    terms: Counter[str] = Counter()
    live = set(wiki_pages())
    for number, row in rows:
        if len(row) != 3:
            fail(failures, GLOSSARY, "glossary row has wrong field count", number)
            continue
        terms[row[0].casefold()] += 1
        owner_links = list(LINK_RE.finditer(row[2]))
        if len(owner_links) != 1:
            fail(failures, GLOSSARY, "glossary Owner must contain exactly one link", number)
            continue
        target, _ = resolve_link(GLOSSARY, owner_links[0].group(2))
        if target not in live:
            fail(failures, GLOSSARY, "glossary Owner is not a live wiki page", number)
    for term, count in sorted(terms.items()):
        if count > 1:
            fail(failures, GLOSSARY, f"duplicate glossary term: {term}")


def affected_pages(cell: str) -> set[str]:
    """Return affected wiki page paths from contradiction cell."""
    result = set()
    for match in LINK_RE.finditer(cell):
        target, _ = resolve_link(CONTRADICTIONS, match.group(2))
        if target is not None and inside(target, WIKI):
            result.add(rel(target))
    return result


def check_contradictions(failures: list[str]) -> None:
    """Check contradiction rows and reciprocal page markers."""
    header, rows = table(CONTRADICTIONS, "ID")
    if header != CONTRADICTION_COLUMNS:
        fail(failures, CONTRADICTIONS, "contradiction table columns do not match schema")
        return
    active: dict[str, tuple[str, set[str]]] = {}
    identifiers: Counter[str] = Counter()
    for number, row in rows:
        if len(row) != len(header):
            fail(failures, CONTRADICTIONS, "contradiction row has wrong field count", number)
            continue
        identifier = row[0].strip("` ")
        identifiers[identifier] += 1
        if not re.fullmatch(r"CONTR-\d{4}", identifier):
            fail(failures, CONTRADICTIONS, f"malformed contradiction ID: {identifier!r}", number)
            continue
        status = row[2].lower()
        if status not in {"contested", "stale", "resolved"}:
            fail(failures, CONTRADICTIONS, f"invalid Claim status: {status!r}", number)
        for index, name in enumerate(CONTRADICTION_COLUMNS):
            if name != "Resolved" and not row[index]:
                fail(failures, CONTRADICTIONS, f"empty contradiction field: {name}", number)
        if not valid_date(row[10]):
            fail(failures, CONTRADICTIONS, f"invalid Opened date: {row[10]!r}", number)
        if status == "resolved":
            if not valid_date(row[11]):
                fail(failures, CONTRADICTIONS, "resolved contradiction lacks date", number)
        elif row[11] not in {"-", "n/a"}:
            fail(failures, CONTRADICTIONS, "active contradiction Resolved must be '-' or 'n/a'", number)
        affected = affected_pages(row[6])
        if not affected:
            fail(failures, CONTRADICTIONS, "Affected pages has no wiki links", number)
        if status in {"contested", "stale"}:
            active[identifier] = (status, affected)
    for identifier, count in identifiers.items():
        if count > 1:
            fail(failures, CONTRADICTIONS, f"duplicate contradiction ID: {identifier}")
    marker_re = re.compile(
        r"Claim status:\s*(contested|stale);\s*contradiction:\s*(CONTR-\d{4})[.]",
        re.IGNORECASE,
    )
    found: dict[str, set[str]] = defaultdict(set)
    for path in wiki_pages():
        if in_isolated(path):
            continue
        if path.name == "index.md" or path in {item.resolve() for item in SUPPORT_PAGES}:
            continue
        text = mask_code(read(path))
        for marker in marker_re.finditer(text):
            status, identifier = marker.group(1).lower(), marker.group(2)
            found[identifier].add(rel(path))
            if identifier not in active:
                fail(failures, path, f"marker references missing/resolved contradiction: {identifier}")
                continue
            if status != active[identifier][0]:
                fail(failures, path, f"marker status disagrees with {identifier}")
            backlink = any(
                target == CONTRADICTIONS.resolve() and slug(anchor) == identifier.lower()
                for _, _, _, target, anchor in links(path)
            )
            if not backlink:
                fail(failures, path, f"marker lacks backlink for {identifier}")
    for identifier, (_, expected) in active.items():
        for missing in sorted(expected - found.get(identifier, set())):
            fail(failures, ROOT / missing, f"missing active marker for {identifier}")
        for extra in sorted(found.get(identifier, set()) - expected):
            fail(failures, ROOT / extra, f"marker page absent from Affected pages for {identifier}")


# Additive GRHayLHD profile. Keep constants and parsing isolated so Illinois
# checks above retain their established behavior.
GR_ROOT = ROOT / "GRHayLHD"
GR_WIKI = WIKI / "grhaylhd"
GR_REGISTRY = ROOT / "raw/grhaylhd/SOURCES.md"
GR_CATALOG = GR_WIKI / "catalog.md"
GR_GLOSSARY = GR_WIKI / "glossary.md"
GR_SOURCE_MAP = GR_WIKI / "source-map.md"
GR_ISSUES = GR_WIKI / "contradictions.md"
GR_SCHEMA = GR_WIKI / "SCHEMA.md"
GR_FOUNDATION = {
    "SCHEMA.md",
    "lint/CHECKS.md",
    "catalog.md",
    "source-map.md",
}
GR_TARGET = GR_FOUNDATION | {
    "index.md",
    "glossary.md",
    "workflows.md",
    "contradictions.md",
    "architecture/index.md",
    "architecture/purpose-build-surface.md",
    "architecture/variables-and-storage.md",
    "architecture/schedule-lifecycle.md",
    "evolution/index.md",
    "evolution/eos-entropy-variants.md",
    "evolution/primitive-conservative-conversion.md",
    "evolution/conservative-recovery.md",
    "evolution/rhs-fluxes-and-sources.md",
    "evolution/matter-boundaries-and-symmetry.md",
    "evolution/perturbations-and-diagnostics.md",
    "integration/index.md",
    "integration/hydrobase-velocity-conversion.md",
    "integration/grhaylib-contract.md",
    "integration/adm-mol-tmunu-contracts.md",
    "integration/parameters-and-configurations.md",
    "validation/index.md",
    "validation/test-inventory-and-oracles.md",
    "validation/coverage-gaps.md",
}
GR_ROOT_ROUTES = {
    "index.md",
    "architecture/index.md",
    "evolution/index.md",
    "integration/index.md",
    "validation/index.md",
    "SCHEMA.md",
    "catalog.md",
    "glossary.md",
    "workflows.md",
    "source-map.md",
    "contradictions.md",
    "lint/CHECKS.md",
}
GR_LEAF_SECTIONS = [
    "Scope and Non-Scope",
    "Summary",
    "Variant Applicability",
    "Claim-Evidence",
    "Details",
    "Caveats",
    "Sources",
    "Related Pages",
]
GR_PAGE_STATUSES = {"draft", "reviewed", "router"}
GR_CLAIM_STATUSES = {
    "declared",
    "visible-implementation",
    "checked-in-observation",
    "unresolved",
    "coverage-gap",
    "out-of-scope",
}
GR_CLAIM_KINDS = {
    "stated-purpose",
    "cactus-interface",
    "parameter",
    "schedule-intent",
    "build-surface",
    "visible-dataflow",
    "visible-formula",
    "visible-call-order",
    "shipped-configuration",
    "test-declaration",
    "numeric-observation",
    "external-behavior",
}
GR_APPLICABILITY = {
    "Common",
    "Hybrid/Simple",
    "Hybrid/Simple+Entropy",
    "Tabulated",
    "Tabulated+Entropy",
}
GR_ISSUE_KINDS = {
    "contradiction",
    "mismatch",
    "hazard",
    "lifecycle-ambiguity",
    "provenance-ambiguity",
}
GR_ISSUE_STATUSES = {"open", "accepted", "resolved"}
GR_CATALOG_COLUMNS = [
    "Page ID",
    "Page",
    "Type",
    "Status",
    "Review date",
    "Aliases / query terms",
]
GR_REGISTRY_COLUMNS = [
    "Source ID",
    "Provenance",
    "Provenance class",
    "Lifecycle",
    "Ingest",
    "Family",
]
GR_PROVENANCE_CLASSES = {
    "narrative",
    "declaration",
    "build input",
    "implementation",
    "shipped configuration",
    "authored test input",
    "checked-in companion configuration",
    "test declaration",
    "checked-in numeric observation",
}
GR_EDGE_COLUMNS = [
    "Source ID",
    "Page ID",
    "Claim kind",
    "Typed locator",
    "Claim status",
    "Next action",
]
GR_ISSUE_COLUMNS = [
    "ID",
    "Kind",
    "Status",
    "Claim / ambiguity",
    "Locator A",
    "Locator B",
    "Affected Page IDs",
    "Safe wording",
    "Impact",
    "Owner / trigger",
    "Resolution test",
    "Resolution locator",
    "Opened",
    "Resolved",
]

# Active-profile knobs used by the shared GRHayL* checker. Defaults preserve
# the established GRHayLHD diagnostics and grammar byte-for-byte.
GR_CHECK_PREFIX = "GRH"
GR_PAGE_PREFIX = "grhaylhd"
GR_TREE_NAME = "GRHayLHD"
GR_SOURCE_ID_PATTERN = r"grhaylhd-[a-z0-9-]+"
GR_ISSUE_ID_PATTERN = r"GRH-\d{4}"
GR_LOCATOR_TYPES = ("doc", "ccl", "build", "c", "macro", "par", "test", "oracle")
GR_LOCATOR_ADMITTED = {
    "doc": {"section"},
    "ccl": {"group", "parameter", "schedule", "storage", "implementation", "requirement"},
    "build": {"field"},
    "c": {"symbol", "call"},
    "macro": {"name", "include"},
    "par": {"parameter", "file"},
    "test": {"case"},
    "oracle": {"dataset", "file"},
}
GR_LOCATOR_QUERIES = {
    ("doc", "section"): {"occurrence"},
    ("ccl", "group"): {"occurrence"},
    ("ccl", "parameter"): {"occurrence"},
    ("ccl", "schedule"): {"context", "occurrence"},
    ("ccl", "storage"): {"occurrence"},
    ("ccl", "implementation"): {"occurrence"},
    ("ccl", "requirement"): {"occurrence"},
    ("build", "field"): {"occurrence"},
    ("c", "symbol"): set(),
    ("c", "call"): {"function"},
    ("macro", "name"): set(),
    ("macro", "include"): set(),
    ("par", "parameter"): {"occurrence"},
    ("par", "file"): set(),
    ("test", "case"): {"occurrence"},
    ("oracle", "dataset"): set(),
    ("oracle", "file"): set(),
}
GR_PARAMETER_REFERENCE_PATTERN = r"USES\s+(?:KEYWORD|CCTK_\w+)"
GR_FORBID_MAGNETIC_NAMES = True
GR_ROOT_INDEX_CHILDREN = {
    "architecture/index.md", "evolution/index.md", "integration/index.md",
    "validation/index.md", "SCHEMA.md", "catalog.md", "glossary.md",
    "workflows.md", "source-map.md", "contradictions.md", "lint/CHECKS.md",
}
GR_APPLICABILITY_LABEL = "variant"

GID_ROOT = ROOT / "GRHayLID"
GID_WIKI = WIKI / "grhaylid"
GID_REGISTRY = ROOT / "raw/grhaylid/SOURCES.md"
GID_FOUNDATION = {"SCHEMA.md", "lint/CHECKS.md", "catalog.md", "source-map.md"}
GID_TARGET = GID_FOUNDATION | {
    "index.md", "glossary.md", "workflows.md", "contradictions.md",
    "architecture/index.md", "architecture/purpose-and-build-surface.md",
    "architecture/schedule-lifecycle.md", "initial-data/index.md",
    "initial-data/one-d-tests-hydro.md", "initial-data/one-d-tests-magnetic.md",
    "initial-data/isotropic-gas.md", "initial-data/constant-density-sphere.md",
    "initial-data/beta-equilibrium.md", "initial-data/entropy-computation.md",
    "integration/index.md", "integration/hydrobase-keyword-extensions.md",
    "integration/grhaylib-contract.md", "integration/parameters-and-configurations.md",
    "validation/index.md", "validation/coverage-gaps.md",
}
GID_ROOT_ROUTES = {
    "index.md", "architecture/index.md", "initial-data/index.md",
    "integration/index.md", "validation/index.md", "SCHEMA.md", "catalog.md",
    "glossary.md", "workflows.md", "source-map.md", "contradictions.md",
    "lint/CHECKS.md",
}
GID_LEAF_SECTIONS = [
    "Scope and Non-Scope", "Summary", "Mode Applicability", "Claim-Evidence",
    "Details", "Caveats", "Sources", "Related Pages",
]
GID_CLAIM_STATUSES = {
    "declared", "visible-implementation", "unresolved", "coverage-gap",
    "out-of-scope",
}
GID_CLAIM_KINDS = {
    "stated-purpose", "cactus-interface", "parameter", "schedule-intent",
    "build-surface", "visible-dataflow", "visible-formula",
    "visible-call-order", "external-behavior",
}
GID_APPLICABILITY = {
    "Common", "HydroTest1D", "HydroTest1D+Magnetic", "IsotropicGas",
    "ConstantDensitySphere", "BetaEquilibrium", "Entropy/Hybrid",
    "Entropy/Tabulated",
}
GID_PROVENANCE_CLASSES = {"narrative", "declaration", "build input", "implementation"}
GID_ROOT_INDEX_CHILDREN = {
    "architecture/index.md", "initial-data/index.md", "integration/index.md",
    "validation/index.md", "SCHEMA.md", "catalog.md", "glossary.md",
    "workflows.md", "source-map.md", "contradictions.md", "lint/CHECKS.md",
}
GID_LOCATOR_ADMITTED = {
    "doc": {"section"},
    "ccl": {"parameter", "schedule", "implementation", "requirement"},
    "build": {"field"},
    "c": {"symbol", "call"},
    "macro": {"name", "include"},
}
GID_LOCATOR_QUERIES = {
    key: value for key, value in GR_LOCATOR_QUERIES.items()
    if key[0] in GID_LOCATOR_ADMITTED and key[1] in GID_LOCATOR_ADMITTED[key[0]]
}
GID_SIBLING_PATH_RES = tuple(
    re.compile(
        r"(?<![\w/.])(?:(?:\.\./)+|\./)?"
        + rf"(?P<source>{re.escape(tree)}/[A-Za-z0-9_./*?+\[\]-]+)"
    )
    for tree in (
        "IllinoisGRMHD", "GRHayLHD", "GRHayLHDX", "GRHayLIDX",
        "NRPyLeakageET",
    )
)


def gr_fail(
    failures: list[str], check_id: str, path: Path, message: str, line: int | None = None
) -> None:
    """Append GRHayLHD diagnostic sortable by check ID then path."""
    check_id = GR_CHECK_PREFIX + check_id[3:]
    location = rel(path)
    if line is not None:
        location += f":{line}"
    failures.append(f"[{check_id}] {location}: {message}")


def gr_live_pages() -> dict[str, Path]:
    """Return namespaced Markdown pages by namespace-relative path."""
    if not GR_WIKI.is_dir():
        return {}
    return {
        path.relative_to(GR_WIKI).as_posix(): path.resolve()
        for path in sorted(GR_WIKI.rglob("*.md"))
    }


def gr_expected_page_id(relative_name: str) -> str:
    """Derive canonical Page ID from namespace-relative path."""
    special_ids = {
        "SCHEMA.md": f"{GR_PAGE_PREFIX}.schema",
        "lint/CHECKS.md": f"{GR_PAGE_PREFIX}.lint",
        "index.md": f"{GR_PAGE_PREFIX}.index",
    }
    if relative_name in special_ids:
        return special_ids[relative_name]
    if relative_name.endswith("/index.md"):
        relative_name = relative_name.removesuffix("/index.md")
    else:
        relative_name = relative_name.removesuffix(".md")
    return GR_PAGE_PREFIX + "." + relative_name.replace("/", ".")


def gr_metadata(path: Path) -> tuple[str | None, str | None]:
    """Extract namespaced page status and review date."""
    head = "\n".join(read(path).splitlines()[:6])
    match = re.search(
        r"^> Page status: ([a-z-]+) · Last reviewed: ([^\s]+)$", head, re.MULTILINE
    )
    return match.groups() if match else (None, None)


def gr_catalog(
    failures: list[str], live: dict[str, Path]
) -> tuple[dict[str, Path], dict[Path, str]]:
    """Validate catalog exact live-page set and metadata."""
    page_ids: dict[str, Path] = {}
    path_ids: dict[Path, str] = {}
    if not GR_CATALOG.is_file():
        gr_fail(failures, "GRH020", GR_CATALOG, "catalog missing")
        return page_ids, path_ids
    header, rows = table(GR_CATALOG, "Page ID")
    if header != GR_CATALOG_COLUMNS:
        gr_fail(failures, "GRH020", GR_CATALOG, "catalog columns do not match schema")
        return page_ids, path_ids
    for number, row in rows:
        if len(row) != len(header):
            gr_fail(failures, "GRH021", GR_CATALOG, "catalog row has wrong field count", number)
            continue
        identifier_match = re.fullmatch(
            rf"`({re.escape(GR_PAGE_PREFIX)}\.[a-z0-9.-]+)`", row[0]
        )
        link_matches = list(LINK_RE.finditer(row[1]))
        if identifier_match is None or len(link_matches) != 1:
            gr_fail(failures, "GRH021", GR_CATALOG, "catalog row needs one valid Page ID and link", number)
            continue
        identifier = identifier_match.group(1)
        target, _ = resolve_link(GR_CATALOG, link_matches[0].group(2))
        if target is None or not inside(target, GR_WIKI):
            gr_fail(failures, "GRH021", GR_CATALOG, "catalog target outside namespace", number)
            continue
        if not target.is_file():
            gr_fail(
                failures, "GRH023", GR_CATALOG,
                f"catalog target is not a live page: {rel(target)}", number,
            )
            continue
        if identifier in page_ids:
            gr_fail(failures, "GRH021", GR_CATALOG, f"duplicate Page ID: {identifier}", number)
        if target in path_ids:
            gr_fail(failures, "GRH021", GR_CATALOG, f"duplicate page: {rel(target)}", number)
        page_ids[identifier] = target
        path_ids[target] = identifier
        relative_name = target.relative_to(GR_WIKI).as_posix()
        expected_id = gr_expected_page_id(relative_name)
        if identifier != expected_id:
            gr_fail(failures, "GRH021", GR_CATALOG, f"Page ID must be {expected_id}", number)
        status, reviewed = gr_metadata(target) if target.is_file() else (None, None)
        expected_type = "router" if target.name == "index.md" else (
            "governance" if target.relative_to(GR_WIKI).as_posix() in {
                "SCHEMA.md", "lint/CHECKS.md", "catalog.md", "source-map.md",
                "workflows.md", "contradictions.md", "glossary.md"
            } else "leaf"
        )
        if row[2] != expected_type:
            gr_fail(failures, "GRH022", GR_CATALOG, f"type disagrees with {rel(target)}", number)
        if row[3] != status or row[4] != reviewed:
            gr_fail(failures, "GRH022", GR_CATALOG, f"metadata disagrees with {rel(target)}", number)
    live_paths = set(live.values())
    for path in sorted(live_paths - set(path_ids)):
        gr_fail(failures, "GRH023", GR_CATALOG, f"missing live page: {rel(path)}")
    for path in sorted(set(path_ids) - live_paths):
        gr_fail(failures, "GRH023", GR_CATALOG, f"catalog target is not live: {rel(path)}")
    return page_ids, path_ids


def gr_glossary(failures: list[str], live: dict[str, Path]) -> None:
    """Validate one unique namespaced owner for every GRHayLHD term."""
    if not GR_GLOSSARY.is_file():
        return
    header, rows = table(GR_GLOSSARY, "Term")
    if header != ["Term", "Routing meaning", "Owner"]:
        gr_fail(failures, "GRH080", GR_GLOSSARY, "glossary columns do not match schema")
        return
    terms: Counter[str] = Counter()
    live_paths = set(live.values())
    for number, row in rows:
        if len(row) != 3 or any(not cell for cell in row):
            gr_fail(failures, "GRH081", GR_GLOSSARY, "invalid glossary row", number)
            continue
        terms[row[0].casefold()] += 1
        owner_links = list(LINK_RE.finditer(row[2]))
        if len(owner_links) != 1:
            gr_fail(failures, "GRH081", GR_GLOSSARY, "Owner needs exactly one link", number)
            continue
        target, _ = resolve_link(GR_GLOSSARY, owner_links[0].group(2))
        if target not in live_paths or target == GR_GLOSSARY.resolve():
            gr_fail(failures, "GRH081", GR_GLOSSARY, "Owner is not another live namespaced page", number)
    for term, count in sorted(terms.items()):
        if count > 1:
            gr_fail(failures, "GRH082", GR_GLOSSARY, f"duplicate glossary term: {term}")


def gr_registry(
    failures: list[str], required: bool
) -> tuple[dict[str, tuple[str, str]], dict[str, str]]:
    """Validate lossless, non-overlapping GRHayLHD source registration."""
    records: dict[str, tuple[str, str]] = {}
    owners: dict[str, str] = {}
    if not GR_REGISTRY.is_file():
        if required:
            gr_fail(failures, "GRH030", GR_REGISTRY, "source registry missing")
        return records, owners
    header, rows = table(GR_REGISTRY, "Source ID")
    if header != GR_REGISTRY_COLUMNS:
        gr_fail(failures, "GRH030", GR_REGISTRY, "registry columns do not match schema")
        return records, owners
    for number, row in rows:
        if len(row) != len(header):
            gr_fail(failures, "GRH031", GR_REGISTRY, "registry row has wrong field count", number)
            continue
        source_match = re.fullmatch(rf"`({GR_SOURCE_ID_PATTERN})`", row[0])
        path_match = re.fullmatch(r"`([^`]+)`", row[1])
        if source_match is None or path_match is None:
            gr_fail(failures, "GRH031", GR_REGISTRY, "Source ID and Provenance need one literal", number)
            continue
        if any(not cell for cell in row):
            gr_fail(failures, "GRH031", GR_REGISTRY, "registry row contains empty cell", number)
            continue
        identifier, pattern = source_match.group(1), path_match.group(1)
        if identifier in records:
            gr_fail(failures, "GRH031", GR_REGISTRY, f"duplicate Source ID: {identifier}", number)
        if row[4] not in {"registered", "partial", "ingested"}:
            gr_fail(failures, "GRH031", GR_REGISTRY, f"invalid ingest state: {row[4]}", number)
        if row[2] not in GR_PROVENANCE_CLASSES or row[3] not in {"living", "frozen"}:
            gr_fail(failures, "GRH031", GR_REGISTRY, "invalid provenance class or lifecycle", number)
        records[identifier] = (pattern, row[4])
        components = pattern.split("/")
        invalid = (
            not pattern or "\\" in pattern or pattern.startswith("/")
            or not pattern.startswith(f"{GR_TREE_NAME}/")
            or any(part in {"", ".", ".."} for part in components)
            or any("**" in part and part != "**" for part in components)
        )
        if invalid:
            gr_fail(failures, "GRH032", GR_REGISTRY, f"invalid or escaping provenance: {pattern}", number)
            continue
        try:
            matches = sorted(ROOT.glob(pattern), key=lambda item: rel(item))
        except (OSError, RuntimeError, ValueError, NotImplementedError) as error:
            gr_fail(failures, "GRH032", GR_REGISTRY, f"invalid provenance {pattern}: {error}", number)
            continue
        if not matches:
            gr_fail(failures, "GRH033", GR_REGISTRY, f"unmatched provenance: {pattern}", number)
            continue
        nonfiles = [item for item in matches if not item.is_file()]
        if GR_PAGE_PREFIX == "grhaylid" and nonfiles:
            gr_fail(
                failures, "GRH033", GR_REGISTRY,
                f"provenance expansion includes non-file: {pattern}", number,
            )
            continue
        files = [item for item in matches if item.is_file()]
        if not files:
            gr_fail(failures, "GRH033", GR_REGISTRY, f"directory-only provenance: {pattern}", number)
            continue
        for item in files:
            if not inside(item, GR_ROOT):
                gr_fail(failures, "GRH032", GR_REGISTRY, f"provenance escapes {GR_TREE_NAME}: {rel(item)}", number)
                continue
            item_name = rel(item)
            if item_name in owners:
                gr_fail(
                    failures, "GRH034", GR_REGISTRY,
                    f"overlap for {item_name}: {owners[item_name]} and {identifier}", number,
                )
            else:
                owners[item_name] = identifier
    command = subprocess.run(
        ["rg", "--files", GR_TREE_NAME], cwd=ROOT, text=True,
        capture_output=True, check=False,
    )
    if command.returncode not in {0, 1}:
        gr_fail(failures, "GRH035", GR_ROOT, f"rg --files failed: {command.stderr.strip()}")
        actual: set[str] = set()
    else:
        actual = set(filter(None, command.stdout.splitlines()))
    for path in sorted(actual - set(owners)):
        gr_fail(failures, "GRH035", ROOT / path, f"{GR_TREE_NAME} file is unregistered")
    for path in sorted(set(owners) - actual):
        gr_fail(failures, "GRH035", ROOT / path, "registry expansion is outside rg inventory")
    return records, owners


def gr_normalize_heading(value: str) -> str:
    """Apply schema heading normalization."""
    return re.sub(r"[ \t]+", " ", unicodedata.normalize("NFKC", value).strip())


def gr_doc_headings(path: Path) -> list[str]:
    """Extract admitted Markdown, setext, README, and TeX headings."""
    text = mask_code(read(path))
    headings = [match.group(1).strip() for match in re.finditer(r"^#{1,6}\s+(.+?)\s*#*\s*$", text, re.MULTILINE)]
    lines = text.splitlines()
    for index in range(len(lines) - 1):
        if lines[index].strip() and re.fullmatch(r"[=-]{3,}", lines[index + 1].strip()):
            headings.append(lines[index].strip())
    headings.extend(
        match.group(1).strip()
        for match in re.finditer(r"^\s*\\(?:part|chapter|section|subsection|subsubsection)\*?\{([^}]+)\}", text, re.MULTILINE)
    )
    return headings


def gr_ccl_schedule_container(text: str, match: re.Match[str]) -> str | None:
    """Return exact CCL schedule bin/group from one declaration line."""
    line_end = text.find("\n", match.start())
    declaration = text[match.start():line_end if line_end >= 0 else len(text)]
    container = re.search(
        r"\b(?:in|at)\s+([A-Za-z_]\w*)\b", declaration, re.IGNORECASE
    )
    return container.group(1) if container else None


def gr_mask_c_literals(text: str) -> str:
    """Mask C comments and literals while retaining offsets and newlines."""
    pattern = re.compile(
        r"//[^\n]*|/\*.*?(?:\*/|$)|\"(?:\\.|[^\"\\])*\"|"
        r"'(?:\\.|[^'\\])*'",
        re.DOTALL,
    )
    return pattern.sub(
        lambda match: "".join(
            "\n" if char == "\n" else " " for char in match.group(0)
        ),
        text,
    )


def gr_occurrence_count(path: Path, locator_type: str, key: str, value: str, query: dict[str, str]) -> int:
    """Count exact stable-locator targets in one local file."""
    text = read(path)
    escaped = re.escape(value)
    if locator_type == "doc" and key == "section":
        wanted = gr_normalize_heading(value)
        return sum(gr_normalize_heading(item) == wanted for item in gr_doc_headings(path))
    if locator_type == "ccl":
        if key == "group":
            pattern = rf"^\s*CCTK_[A-Z0-9_]+\s+{escaped}\s+type\s*="
        elif key == "parameter":
            pattern = (
                rf"^\s*(?:(?:INT|REAL|KEYWORD|STRING|"
                rf"CCTK_(?:BOOLEAN|INT|REAL|STRING))\s+|"
                rf"{GR_PARAMETER_REFERENCE_PATTERN}\s+){escaped}\b"
            )
        elif key == "schedule":
            pattern = rf"^\s*schedule(?:\s+group)?\s+{escaped}\b"
        elif key == "storage":
            pattern = rf"^\s*STORAGE:\s*[^\n]*\b{escaped}(?:\[\d+\])?(?=\s*(?:,|$))"
        elif key == "implementation":
            pattern = rf"^\s*IMPLEMENTS:\s*{escaped}\s*$"
        elif key == "requirement":
            pattern = rf"^\s*requires\s+{escaped}\s*$"
        else:
            return 0
        matches = list(re.finditer(pattern, text, re.MULTILINE | re.IGNORECASE))
        context = query.get("context")
        if context:
            if key != "schedule":
                return 0
            matches = [
                match
                for match in matches
                if (gr_ccl_schedule_container(text, match) or "").casefold()
                == context.casefold()
            ]
        return len(matches)
    if locator_type == "build" and key == "field":
        return len(re.findall(rf"^\s*{escaped}\s*=", text, re.MULTILINE))
    if locator_type == "macro" and key == "name":
        return len(re.findall(rf"^\s*#\s*define\s+{escaped}\b", text, re.MULTILINE))
    if locator_type == "macro" and key == "include":
        return len(re.findall(rf'^\s*#\s*include\s+[<"]{escaped}[>"]', text, re.MULTILINE))
    if locator_type == "c" and key == "symbol":
        symbol_text = (
            gr_mask_c_literals(text)
            if GR_PAGE_PREFIX == "grhaylid" else text
        )
        pattern = rf"^\s*(?!if\b|for\b|while\b|switch\b)(?:[A-Za-z_]\w*[\s*]+)+{escaped}\s*\([^;]*\)\s*\{{"
        return len(re.findall(pattern, symbol_text, re.MULTILINE))
    if locator_type == "c" and key == "call":
        function = query.get("function")
        if not function:
            return 0
        masked = gr_mask_c_literals(text)
        function_start = re.search(
            rf"^\s*(?:[A-Za-z_]\w*[\s*]+)+{re.escape(function)}"
            rf"\s*\([^;]*\)\s*\{{",
            masked,
            re.MULTILINE,
        )
        if function_start is None:
            return 0
        depth = 1
        body_end = len(masked)
        for index in range(function_start.end(), len(masked)):
            if masked[index] == "{":
                depth += 1
            elif masked[index] == "}":
                depth -= 1
                if depth == 0:
                    body_end = index
                    break
        return len(re.findall(
            rf"\b{escaped}\s*\(", masked[function_start.end():body_end]
        ))
    if locator_type == "par" and key == "parameter":
        return len(re.findall(rf"^\s*{escaped}\s*=", text, re.MULTILINE | re.IGNORECASE))
    if locator_type == "par" and key == "file" and value == "file":
        return 1
    if locator_type == "test" and key == "case":
        return len(re.findall(rf"^\s*TEST\s+{escaped}\b", text, re.MULTILINE | re.IGNORECASE))
    if locator_type == "oracle" and key == "dataset":
        return len(re.findall(rf"(?:^|\s)(?:dataset|name)\s*[:=]\s*{escaped}(?:\s|$)", text, re.MULTILINE | re.IGNORECASE))
    if locator_type == "oracle" and key == "file" and value == "file":
        return 1
    return 0


def gr_resolve_locator(
    failures: list[str], locator: str, owners: dict[str, str],
    expected_source: str | None, source_path: Path, line: int | None = None,
) -> bool:
    """Resolve one typed locator and enforce registry ownership/uniqueness."""
    locator_types = "|".join(re.escape(item) for item in GR_LOCATOR_TYPES)
    match = re.fullmatch(
        rf"({locator_types}):({re.escape(GR_TREE_NAME)}/[^#?]+)"
        rf"#([^?]+)(?:\?(.*))?",
        locator,
    )
    if match is None:
        gr_fail(failures, "GRH040", source_path, f"malformed typed locator: {locator}", line)
        return False
    locator_type, relative, fragment, query_text = match.groups()
    target = (ROOT / relative).resolve()
    if not target.is_file() or not inside(target, GR_ROOT):
        gr_fail(failures, "GRH041", source_path, f"locator path missing/outside scope: {relative}", line)
        return False
    owner = owners.get(relative)
    if owner is None:
        gr_fail(failures, "GRH041", source_path, f"locator path is unregistered: {relative}", line)
        return False
    if expected_source is not None and owner != expected_source:
        gr_fail(failures, "GRH041", source_path, f"locator owned by {owner}, not {expected_source}", line)
        return False
    query: dict[str, str] = {}
    if query_text:
        for item in query_text.split("&"):
            if "=" not in item:
                gr_fail(failures, "GRH040", source_path, f"malformed locator query: {locator}", line)
                return False
            key, value = item.split("=", 1)
            if key in query or key not in {"occurrence", "context", "function"} or not value:
                gr_fail(failures, "GRH040", source_path, f"invalid locator query: {locator}", line)
                return False
            query[key] = value
    if fragment == "file":
        key, value = "file", "file"
    elif "=" in fragment:
        key, value = fragment.split("=", 1)
        if key == "file":
            gr_fail(
                failures, "GRH040", source_path,
                f"file qualifier must be exact '#file': {locator}", line,
            )
            return False
    else:
        gr_fail(failures, "GRH040", source_path, f"malformed locator qualifier: {locator}", line)
        return False
    if key not in GR_LOCATOR_ADMITTED[locator_type] or not value:
        gr_fail(failures, "GRH040", source_path, f"qualifier not admitted for {locator_type}: {fragment}", line)
        return False
    if set(query) - GR_LOCATOR_QUERIES[(locator_type, key)]:
        gr_fail(failures, "GRH040", source_path, f"query not admitted for qualifier: {locator}", line)
        return False
    if locator_type == "ccl" and "context" in query and "occurrence" in query:
        gr_fail(failures, "GRH040", source_path, f"context and occurrence are mutually exclusive: {locator}", line)
        return False
    count = gr_occurrence_count(target, locator_type, key, value, query)
    occurrence = query.get("occurrence")
    if occurrence:
        if not occurrence.isdigit() or int(occurrence) < 1 or int(occurrence) > count:
            gr_fail(failures, "GRH042", source_path, f"locator occurrence does not resolve: {locator}", line)
            return False
        count = 1
    if count != 1:
        adjective = "unresolved" if count == 0 else "non-unique"
        gr_fail(failures, "GRH042", source_path, f"{adjective} typed locator: {locator}", line)
        return False
    return True


def gr_pages(failures: list[str], live: dict[str, Path]) -> None:
    """Validate namespaced page inventory, metadata, and leaf/router schema."""
    for missing in sorted(GR_FOUNDATION - set(live)):
        gr_fail(failures, "GRH010", GR_WIKI / missing, "required foundation page missing")
    for extra in sorted(set(live) - GR_TARGET):
        gr_fail(failures, "GRH010", GR_WIKI / extra, "page outside exact target tree")
    for name, path in sorted(live.items()):
        if GR_FORBID_MAGNETIC_NAMES and "magnetic" in name.casefold():
            gr_fail(failures, "GRH011", path, "Magnetics branch forbidden for GRHayLHD")
        status, reviewed = gr_metadata(path)
        if status not in GR_PAGE_STATUSES:
            gr_fail(failures, "GRH012", path, f"invalid or missing page status: {status}")
        if status == "router":
            if reviewed != "n/a":
                gr_fail(failures, "GRH012", path, "router review date must be n/a")
        elif not valid_date(reviewed or ""):
            gr_fail(failures, "GRH012", path, f"invalid review date: {reviewed}")
        is_domain_leaf = name not in {
            "SCHEMA.md", "lint/CHECKS.md", "catalog.md", "source-map.md",
            "workflows.md", "contradictions.md", "glossary.md",
        } and not name.endswith("/index.md") and name != "index.md"
        raw_text = read(path)
        text = mask_code(raw_text)
        if is_domain_leaf:
            raw_lines = raw_text.splitlines()
            masked_lines = text.splitlines()
            atx_h1_count = len(re.findall(
                r"^[ ]{0,3}#(?!#)(?:[ \t]+.*)?$", text, re.MULTILINE
            ))
            setext_h1_count = sum(
                bool(masked_lines[index].strip())
                and bool(re.fullmatch(
                    r"[ ]{0,3}=+[ \t]*", masked_lines[index + 1]
                ))
                for index in range(len(masked_lines) - 1)
            )
            h1_count = atx_h1_count + setext_h1_count
            up_match = (
                re.fullmatch(r"> Up: \[[^]]+\]\(([^)]+)\)", raw_lines[3])
                if len(raw_lines) > 3 else None
            )
            up_is_valid = up_match is not None
            if GR_PAGE_PREFIX == "grhaylid" and up_match is not None:
                up_is_valid = (
                    resolve_link(path, up_match.group(1))[0]
                    == (path.parent / "index.md").resolve()
                )
            exact_leaf_preamble = (
                len(raw_lines) >= 5
                and h1_count == 1
                and bool(re.fullmatch(r"#(?!#)\s+.+", raw_lines[0]))
                and raw_lines[1] == ""
                and bool(re.fullmatch(
                    r"> Page status: (?:draft|reviewed) · Last reviewed: "
                    r"\d{2}-\d{2}-\d{4}",
                    raw_lines[2],
                ))
                and up_is_valid
                and raw_lines[4] == ""
            )
            if not exact_leaf_preamble:
                gr_fail(
                    failures, "GRH013", path,
                    "leaf preamble disagrees with schema",
                )
            headings = [match.group(1).strip() for match in re.finditer(r"^##\s+(.+?)\s*$", text, re.MULTILINE)]
            if headings != GR_LEAF_SECTIONS:
                gr_fail(failures, "GRH013", path, "leaf H2 sections disagree with schema")
            for heading in GR_LEAF_SECTIONS:
                if not section(text, heading):
                    gr_fail(failures, "GRH013", path, f"leaf section is empty: {heading}")
            if status not in {"draft", "reviewed"}:
                gr_fail(failures, "GRH013", path, "leaf status must be draft or reviewed")
            claim_header, claim_rows = table(path, "Claim ID")
            if claim_header != ["Claim ID", "Claim", "Status", "Evidence", "Typed locator"]:
                gr_fail(failures, "GRH014", path, "claim-evidence table columns do not match schema")
            claim_ids: Counter[str] = Counter()
            for number, row in claim_rows:
                identifier = row[0].strip("` ") if row else ""
                claim_ids[identifier] += 1
                if (
                    len(row) != 5 or any(not cell for cell in row)
                    or not re.fullmatch(r"[A-Z][A-Z0-9-]*-\d+", identifier)
                    or row[2] not in GR_CLAIM_STATUSES
                ):
                    gr_fail(failures, "GRH014", path, "invalid claim-evidence row", number)
            for identifier, count in sorted(claim_ids.items()):
                if count > 1:
                    gr_fail(failures, "GRH014", path, f"duplicate Claim ID: {identifier}")
            applicability_header, applicability_rows = table(path, "Applicability")
            if not applicability_header or not applicability_rows:
                gr_fail(
                    failures, "GRH014", path,
                    f"{GR_APPLICABILITY_LABEL}-applicability table missing or empty",
                )
            for number, row in applicability_rows:
                if not row or row[0] not in GR_APPLICABILITY:
                    gr_fail(
                        failures, "GRH014", path,
                        f"invalid {GR_APPLICABILITY_LABEL} applicability", number,
                    )
        elif name.endswith("index.md") or name == "index.md":
            if GR_PAGE_PREFIX == "grhaylid":
                router_h1s = re.findall(
                    r"^[ ]{0,3}#(?!#)[ \t]+\S.*$", text, re.MULTILINE
                )
                if len(router_h1s) != 1:
                    gr_fail(
                        failures, "GRH015", path,
                        "router must contain exactly one nonempty H1",
                    )
            if status != "router":
                gr_fail(failures, "GRH015", path, "router must use router status")
            up_matches = [
                match
                for line in read(path).splitlines()
                for match in [re.fullmatch(r"> Up: \[[^]]+\]\(([^)]+)\)", line)]
                if match is not None
            ]
            up_target = resolve_link(path, up_matches[0].group(1))[0] if len(up_matches) == 1 else None
            valid_up = (
                up_target == AGENTS.resolve()
                if name == "index.md"
                else up_target is not None and inside(up_target, GR_WIKI)
            )
            if len(up_matches) != 1 or not valid_up:
                gr_fail(failures, "GRH015", path, "router lacks exact namespaced Up target")
            header, route_rows = table(path, "Page")
            if header != ["Page", "Use it for"]:
                gr_fail(failures, "GRH015", path, "router lacks Page / Use it for table")
            expected_children: set[Path]
            if name == "index.md":
                expected_children = {
                    child for child_name, child in live.items()
                    if child_name in GR_ROOT_INDEX_CHILDREN
                }
            else:
                branch = name.removesuffix("/index.md") + "/"
                expected_children = {
                    child for child_name, child in live.items()
                    if child_name.startswith(branch) and child_name != name
                    and "/" not in child_name[len(branch):]
                }
            if GR_PAGE_PREFIX == "grhaylid":
                routed_counts: Counter[Path] = Counter()
                live_paths = set(live.values())
                for number, row in route_rows:
                    matches = list(LINK_RE.finditer(row[0])) if len(row) == 2 else []
                    if len(row) != 2 or any(not cell for cell in row) or len(matches) != 1:
                        gr_fail(failures, "GRH016", path, "invalid router row", number)
                        continue
                    target, _ = resolve_link(path, matches[0].group(2))
                    if target is None or not target.is_file() or target not in live_paths:
                        gr_fail(
                            failures, "GRH016", path,
                            "invalid or missing router target", number,
                        )
                        continue
                    routed_counts[target] += 1
                routed = set(routed_counts)
            else:
                routed = {
                    target
                    for _, row in route_rows
                    if row
                    for match in LINK_RE.finditer(row[0])
                    for target, _ in [resolve_link(path, match.group(2))]
                    if target is not None and target in set(live.values())
                }
            for child in sorted(expected_children - routed):
                gr_fail(failures, "GRH016", path, f"router misses live child: {rel(child)}")
            for child in sorted(routed - expected_children):
                gr_fail(failures, "GRH016", path, f"router table has non-child: {rel(child)}")
            if GR_PAGE_PREFIX == "grhaylid":
                for child, count in sorted(routed_counts.items()):
                    if child in expected_children and count != 1:
                        gr_fail(
                            failures, "GRH016", path,
                            f"router repeats live child: {rel(child)}",
                        )


def gr_edges(
    failures: list[str], page_ids: dict[str, Path], path_ids: dict[Path, str],
    records: dict[str, tuple[str, str]], owners: dict[str, str],
) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Validate canonical source-page edge table and leaf claim locators."""
    source_consumers: dict[str, set[str]] = defaultdict(set)
    page_sources: dict[str, set[str]] = defaultdict(set)
    page_locator_statuses: dict[str, set[tuple[str, str]]] = defaultdict(set)
    if not GR_SOURCE_MAP.is_file():
        return source_consumers, page_sources
    header, rows = table(GR_SOURCE_MAP, "Source ID")
    if header != GR_EDGE_COLUMNS:
        gr_fail(failures, "GRH050", GR_SOURCE_MAP, "source-map columns do not match schema")
        return source_consumers, page_sources
    seen: set[tuple[str, str, str, str, str]] = set()
    for number, row in rows:
        if len(row) != len(header):
            gr_fail(failures, "GRH051", GR_SOURCE_MAP, "edge row has wrong field count", number)
            continue
        if any(not cell for cell in row):
            gr_fail(failures, "GRH051", GR_SOURCE_MAP, "edge row contains empty cell", number)
            continue
        source_id = row[0].strip("` ")
        page_id = row[1].strip("` ")
        locator = row[3].strip("` ")
        key = (source_id, page_id, row[2], locator, row[4])
        if key in seen:
            gr_fail(failures, "GRH051", GR_SOURCE_MAP, "duplicate source-page edge", number)
        seen.add(key)
        if source_id not in records:
            gr_fail(failures, "GRH051", GR_SOURCE_MAP, f"unknown Source ID: {source_id}", number)
        if page_id not in page_ids:
            gr_fail(failures, "GRH051", GR_SOURCE_MAP, f"unknown Page ID: {page_id}", number)
        if row[4] not in GR_CLAIM_STATUSES:
            gr_fail(failures, "GRH051", GR_SOURCE_MAP, f"invalid claim status: {row[4]}", number)
        if row[2] not in GR_CLAIM_KINDS:
            gr_fail(failures, "GRH051", GR_SOURCE_MAP, f"invalid claim kind: {row[2]}", number)
        if source_id in records:
            gr_resolve_locator(failures, locator, owners, source_id, GR_SOURCE_MAP, number)
        if source_id in records and page_id in page_ids:
            source_consumers[source_id].add(page_id)
            page_sources[page_id].add(source_id)
            page_locator_statuses[page_id].add((locator, row[4]))
    for source_id, (_, ingest) in sorted(records.items()):
        if ingest == "ingested" and not source_consumers.get(source_id):
            gr_fail(failures, "GRH052", GR_SOURCE_MAP, f"ingested source lacks consumer: {source_id}")
    target_page_ids = {gr_expected_page_id(name) for name in GR_TARGET}
    if target_page_ids.issubset(page_ids):
        for source_id in sorted(set(records) - set(source_consumers)):
            gr_fail(failures, "GRH052", GR_SOURCE_MAP, f"complete target leaves source unconsumed: {source_id}")
    for path, page_id in sorted(path_ids.items(), key=lambda item: item[1]):
        name = path.relative_to(GR_WIKI).as_posix()
        status, _ = gr_metadata(path)
        is_leaf = name not in {
            "SCHEMA.md", "lint/CHECKS.md", "catalog.md", "source-map.md",
            "workflows.md", "contradictions.md", "glossary.md",
        } and not name.endswith("/index.md") and name != "index.md"
        if is_leaf and status == "reviewed" and not page_sources.get(page_id):
            gr_fail(failures, "GRH053", path, "reviewed leaf lacks canonical source edge")
        if is_leaf:
            _, claim_rows = table(path, "Claim ID")
            for number, row in claim_rows:
                if len(row) == 5:
                    locator = row[4].strip("` ")
                    gr_resolve_locator(failures, locator, owners, None, path, number)
                    if (locator, row[2]) not in page_locator_statuses.get(
                        page_id, set()
                    ):
                        gr_fail(
                            failures, "GRH054", path,
                            "claim locator/status lacks exact canonical page edge",
                            number,
                        )
            for number, _, _, target, _ in links(path):
                if target is None or not target.is_file() or not inside(target, GR_ROOT):
                    continue
                source_id = owners.get(rel(target))
                if source_id is not None and source_id not in page_sources.get(page_id, set()):
                    gr_fail(
                        failures, "GRH055", path,
                        f"direct evidence source lacks canonical page edge: {source_id}", number,
                    )
    return source_consumers, page_sources


def gr_resolve_resolution_locator(
    failures: list[str], locator: str, owners: dict[str, str], line: int
) -> bool:
    """Resolve source evidence or namespaced documentation disposition."""
    if not locator.startswith("kb:"):
        return gr_resolve_locator(
            failures, locator, owners, None, GR_ISSUES, line
        )
    match = re.fullmatch(
        rf"kb:(wiki/{re.escape(GR_PAGE_PREFIX)}/[^#]+)#section=(.+)", locator
    )
    if match is None:
        gr_fail(
            failures, "GRH061", GR_ISSUES,
            f"malformed resolution locator: {locator}", line,
        )
        return False
    target = (ROOT / match.group(1)).resolve()
    if not target.is_file() or not inside(target, GR_WIKI):
        gr_fail(
            failures, "GRH061", GR_ISSUES,
            f"resolution disposition missing/outside namespace: {locator}",
            line,
        )
        return False
    wanted = gr_normalize_heading(match.group(2))
    if sum(
        gr_normalize_heading(item) == wanted for item in gr_doc_headings(target)
    ) != 1:
        gr_fail(
            failures, "GRH061", GR_ISSUES,
            f"resolution disposition does not resolve uniquely: {locator}",
            line,
        )
        return False
    return True


def gr_issues(
    failures: list[str], page_ids: dict[str, Path], owners: dict[str, str]
) -> None:
    """Validate typed issues and reciprocal backlinks for live affected pages."""
    if not GR_ISSUES.is_file():
        return
    header, rows = table(GR_ISSUES, "ID")
    if header != GR_ISSUE_COLUMNS:
        gr_fail(failures, "GRH060", GR_ISSUES, "issue columns do not match schema")
        return
    seen: set[str] = set()
    affected_pages: dict[str, set[str]] = {}
    for number, row in rows:
        if len(row) != len(header):
            gr_fail(failures, "GRH061", GR_ISSUES, "issue row has wrong field count", number)
            continue
        if any(not cell for cell in row):
            gr_fail(failures, "GRH061", GR_ISSUES, "issue row contains empty cell", number)
            continue
        identifier = row[0].strip("` ")
        if not re.fullmatch(GR_ISSUE_ID_PATTERN, identifier) or identifier in seen:
            gr_fail(failures, "GRH061", GR_ISSUES, f"malformed/duplicate issue ID: {identifier}", number)
        seen.add(identifier)
        if row[1] not in GR_ISSUE_KINDS or row[2] not in GR_ISSUE_STATUSES:
            gr_fail(failures, "GRH061", GR_ISSUES, "invalid issue kind or status", number)
        if not valid_date(row[12]):
            gr_fail(failures, "GRH061", GR_ISSUES, "invalid issue opened date", number)
        resolution_locator = row[11].strip("` ")
        if row[2] == "resolved":
            if resolution_locator in {"-", "n/a"}:
                gr_fail(
                    failures, "GRH061", GR_ISSUES,
                    "resolved issue lacks resolution locator", number,
                )
            else:
                gr_resolve_resolution_locator(
                    failures, resolution_locator, owners, number
                )
            if not valid_date(row[13]):
                gr_fail(failures, "GRH061", GR_ISSUES, "resolved issue lacks date", number)
        else:
            if resolution_locator not in {"-", "n/a"}:
                gr_fail(
                    failures, "GRH061", GR_ISSUES,
                    "active issue has resolution locator", number,
                )
            if row[13] not in {"-", "n/a"}:
                gr_fail(
                    failures, "GRH061", GR_ISSUES,
                    "active issue resolved date must be '-' or 'n/a'", number,
                )
        for locator in row[4:6]:
            if locator not in {"-", "n/a"}:
                gr_resolve_locator(failures, locator.strip("` "), owners, None, GR_ISSUES, number)
        affected = re.findall(
            rf"`({re.escape(GR_PAGE_PREFIX)}\.[a-z0-9.-]+)`", row[6]
        )
        affected_pages[identifier] = set(affected)
        if not affected:
            gr_fail(failures, "GRH061", GR_ISSUES, "issue lacks affected Page IDs", number)
        for page_id in affected:
            expected_ids = {gr_expected_page_id(name) for name in GR_TARGET}
            if page_id not in expected_ids:
                gr_fail(failures, "GRH061", GR_ISSUES, f"unknown affected Page ID: {page_id}", number)
                continue
            target = page_ids.get(page_id)
            if target is None:
                continue
            anchor = identifier.casefold()
            backlink = any(
                link_target == GR_ISSUES.resolve() and slug(link_anchor) == anchor
                for _, _, _, link_target, link_anchor in links(target)
            )
            if not backlink:
                gr_fail(failures, "GRH062", target, f"missing issue backlink: {identifier}")
    for page_id, target in sorted(page_ids.items()):
        for _, _, _, link_target, link_anchor in links(target):
            if link_target != GR_ISSUES.resolve() or not link_anchor:
                continue
            identifier = link_anchor.upper()
            if identifier in affected_pages and page_id not in affected_pages[identifier]:
                gr_fail(failures, "GRH063", target, f"issue backlink absent from affected pages: {identifier}")


def gr_navigation(failures: list[str], live: dict[str, Path]) -> str:
    """Enforce namespace isolation and atomic publication gate."""
    publication = "unpublished"
    index = GR_WIKI / "index.md"
    if index.is_file():
        match = re.search(r"^> Publication: ([a-z]+)$", read(index), re.MULTILINE)
        if match is None:
            gr_fail(failures, "GRH070", index, "missing publication state")
        else:
            publication = match.group(1)
        if publication not in {"unpublished", "ready"}:
            gr_fail(failures, "GRH070", index, f"invalid publication state: {publication}")
    root_targets = {
        target for _, _, _, target, _ in links(AGENTS)
        if target is not None and (
            inside(target, GR_WIKI)
            or (GR_PAGE_PREFIX == "grhaylid" and target == GR_REGISTRY.resolve())
        )
    }
    if publication != "ready" and root_targets:
        gr_fail(failures, "GRH071", AGENTS, f"root {GR_TREE_NAME} route exists while branch is unpublished")
    if publication == "ready":
        if GR_PAGE_PREFIX == "grhaylid":
            expected_router_targets = {
                (GR_WIKI / item).resolve()
                for item in GR_ROOT_ROUTES - {"index.md"}
            } | {GR_REGISTRY.resolve()}
            router_valid, router_rows = named_section_table(
                AGENTS, "GRHayLID Router", ["Go to", "Use it for"]
            )
            if not router_valid:
                gr_fail(
                    failures, "GRH071", AGENTS,
                    "ready branch lacks exact GRHayLID Router table",
                )
            router_targets: Counter[Path] = Counter()
            for number, row in router_rows:
                match = exact_link_match(row[0]) if len(row) == 2 else None
                if len(row) != 2 or any(not cell for cell in row) or match is None:
                    gr_fail(
                        failures, "GRH071", AGENTS,
                        "invalid GRHayLID Router row", number,
                    )
                    continue
                target, _ = resolve_link(AGENTS, match.group(2))
                if target is None or target not in expected_router_targets:
                    gr_fail(
                        failures, "GRH071", AGENTS,
                        "invalid GRHayLID Router target", number,
                    )
                    continue
                router_targets[target] += 1
            for target in sorted(expected_router_targets - set(router_targets)):
                gr_fail(
                    failures, "GRH071", AGENTS,
                    f"ready branch lacks root route: {rel(target)}",
                )
            for target, count in sorted(router_targets.items()):
                if count != 1:
                    gr_fail(
                        failures, "GRH071", AGENTS,
                        f"ready branch repeats root route: {rel(target)}",
                    )

            selector_valid, selector_rows = named_section_table(
                AGENTS,
                "Branch Selector",
                ["Branch", "Exclusive evidence scope", "Start here"],
            )
            gid_selector_rows = [
                (number, row)
                for number, row in selector_rows
                if row and "GRHayLID" in row[0]
            ]
            if not selector_valid or len(gid_selector_rows) != 1:
                gr_fail(
                    failures, "GRH071", AGENTS,
                    "ready branch needs exactly one GRHayLID Branch Selector row",
                )
            else:
                number, row = gid_selector_rows[0]
                match = exact_link_match(row[2]) if len(row) == 3 else None
                target = (
                    resolve_link(AGENTS, match.group(2))[0]
                    if match is not None else None
                )
                if (
                    len(row) != 3
                    or row[0] != "GRHayLID"
                    or row[1] != "`GRHayLID/**`"
                    or match is None
                    or target != index.resolve()
                ):
                    gr_fail(
                        failures, "GRH071", AGENTS,
                        "invalid GRHayLID Branch Selector row", number,
                    )
        elif index.resolve() not in root_targets:
            gr_fail(failures, "GRH071", AGENTS, "ready branch lacks root index route")
    allowed_governance = {
        AGENTS.resolve(), GR_REGISTRY.resolve(), Path(__file__).resolve(),
        (WIKI / "lint/CHECKS.md").resolve(),
    }
    for name, path in sorted(live.items()):
        if GR_PAGE_PREFIX == "grhaylhd":
            forbidden_mentions = source_path_mentions(read(path))
        else:
            masked = mask_code(read(path))
            forbidden_mentions = [
                (
                    masked.count("\n", 0, match.start()) + 1,
                    match.group("source").rstrip(".,;:"),
                )
                for pattern in GID_SIBLING_PATH_RES
                for match in pattern.finditer(masked)
            ]
        for number, mentioned in forbidden_mentions:
            message = (
                "IllinoisGRMHD source citation forbidden in GRHayLHD branch: "
                if GR_PAGE_PREFIX == "grhaylhd"
                else f"sibling source citation forbidden in {GR_TREE_NAME} branch: "
            )
            gr_fail(
                failures, "GRH072", path,
                f"{message}{mentioned}",
                number,
            )
        for number, _, raw, target, _ in links(path):
            if target is None or inside(target, GR_WIKI) or inside(target, GR_ROOT):
                continue
            if target in allowed_governance and (target != AGENTS.resolve() or name == "index.md"):
                continue
            gr_fail(failures, "GRH072", path, f"navigation/citation escapes namespace: {raw}", number)
    if publication == "ready":
        missing = GR_TARGET - set(live)
        for name in sorted(missing):
            gr_fail(failures, "GRH073", GR_WIKI / name, "ready branch target missing")
        for name, path in sorted(live.items()):
            status, _ = gr_metadata(path)
            if name.endswith("/index.md") or name == "index.md":
                if status != "router":
                    gr_fail(failures, "GRH073", path, "ready router status invalid")
            elif status != "reviewed":
                gr_fail(failures, "GRH073", path, "ready branch contains unreviewed page")
        graph: dict[Path, set[Path]] = defaultdict(set)
        for path in live.values():
            for _, _, _, target, _ in links(path):
                if target in set(live.values()):
                    graph[path].add(target)
        distance: dict[Path, int] = {index.resolve(): 0}
        queue = deque([index.resolve()])
        while queue:
            current = queue.popleft()
            for target in sorted(graph.get(current, set())):
                if target not in distance:
                    distance[target] = distance[current] + 1
                    queue.append(target)
        for name, path in sorted(live.items()):
            if path.resolve() not in distance:
                gr_fail(failures, "GRH074", path, "ready page unreachable from branch index")
            elif name not in {"SCHEMA.md", "lint/CHECKS.md", "catalog.md", "source-map.md", "workflows.md", "contradictions.md", "glossary.md"} and distance[path.resolve()] > 2:
                gr_fail(failures, "GRH074", path, "domain page is more than two hops from branch index")
    return publication


def check_grhaylhd(failures: list[str]) -> None:
    """Run isolated additive GRHayLHD profile when namespace exists."""
    if not GR_WIKI.exists() and not GR_REGISTRY.exists():
        return
    live = gr_live_pages()
    gr_pages(failures, live)
    page_ids, path_ids = gr_catalog(failures, live)
    gr_glossary(failures, live)
    required_registry = GR_REGISTRY.exists() or any(
        name in live for name in {"workflows.md", "contradictions.md"}
    )
    records, owners = gr_registry(failures, required_registry)
    gr_edges(failures, page_ids, path_ids, records, owners)
    gr_issues(failures, page_ids, owners)
    gr_navigation(failures, live)


def check_grhaylid(failures: list[str]) -> None:
    """Run isolated additive GRHayLID profile when namespace exists."""
    if not GID_WIKI.exists() and not GID_REGISTRY.exists():
        return
    profile = {
        "GR_ROOT": GID_ROOT,
        "GR_WIKI": GID_WIKI,
        "GR_REGISTRY": GID_REGISTRY,
        "GR_CATALOG": GID_WIKI / "catalog.md",
        "GR_GLOSSARY": GID_WIKI / "glossary.md",
        "GR_SOURCE_MAP": GID_WIKI / "source-map.md",
        "GR_ISSUES": GID_WIKI / "contradictions.md",
        "GR_SCHEMA": GID_WIKI / "SCHEMA.md",
        "GR_FOUNDATION": GID_FOUNDATION,
        "GR_TARGET": GID_TARGET,
        "GR_ROOT_ROUTES": GID_ROOT_ROUTES,
        "GR_LEAF_SECTIONS": GID_LEAF_SECTIONS,
        "GR_CLAIM_STATUSES": GID_CLAIM_STATUSES,
        "GR_CLAIM_KINDS": GID_CLAIM_KINDS,
        "GR_APPLICABILITY": GID_APPLICABILITY,
        "GR_PROVENANCE_CLASSES": GID_PROVENANCE_CLASSES,
        "GR_CHECK_PREFIX": "GID",
        "GR_PAGE_PREFIX": "grhaylid",
        "GR_TREE_NAME": "GRHayLID",
        "GR_SOURCE_ID_PATTERN": r"grhaylid-[a-z0-9-]+",
        "GR_ISSUE_ID_PATTERN": r"GID-\d{4}",
        "GR_LOCATOR_TYPES": ("doc", "ccl", "build", "c", "macro"),
        "GR_LOCATOR_ADMITTED": GID_LOCATOR_ADMITTED,
        "GR_LOCATOR_QUERIES": GID_LOCATOR_QUERIES,
        "GR_PARAMETER_REFERENCE_PATTERN": (
            r"(?:USES|EXTENDS)\s+(?:KEYWORD|CCTK_\w+)"
        ),
        "GR_FORBID_MAGNETIC_NAMES": False,
        "GR_ROOT_INDEX_CHILDREN": GID_ROOT_INDEX_CHILDREN,
        "GR_APPLICABILITY_LABEL": "mode",
    }
    previous = {name: globals()[name] for name in profile}
    globals().update(profile)
    try:
        live = gr_live_pages()
        gr_pages(failures, live)
        page_ids, path_ids = gr_catalog(failures, live)
        gr_glossary(failures, live)
        required_registry = GR_REGISTRY.exists() or any(
            name in live for name in {"workflows.md", "contradictions.md"}
        )
        records, owners = gr_registry(failures, required_registry)
        gr_edges(failures, page_ids, path_ids, records, owners)
        gr_issues(failures, page_ids, owners)
        gr_navigation(failures, live)
    finally:
        globals().update(previous)


def main() -> int:
    """Run one deterministic dispatcher for default and --all modes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all",
        action="store_true",
        help="compatibility alias; identical deterministic coverage",
    )
    parser.parse_args()

    failures: list[str] = []
    check_inventory(failures)
    required = {AGENTS, SOURCES, Path(__file__).resolve()} | {
        WIKI / name for name in EXPECTED_WIKI
    }
    if any(not path.is_file() for path in required):
        print("KB lint failed:")
        for item in sorted(set(failures)):
            print(f"- {item}")
        return 1

    check_text(failures)
    graph = check_links(failures)
    check_routers(failures)
    check_leaves(failures)
    identifiers, covered = parse_manifest(failures)
    check_source_citations(failures, identifiers, covered)
    check_reachability(failures, graph)
    check_catalog(failures)
    check_source_map(failures, identifiers, covered)
    check_glossary(failures)
    check_contradictions(failures)
    check_grhaylhd(failures)
    check_grhaylid(failures)

    if failures:
        print("KB lint failed:")
        for item in sorted(set(failures)):
            print(f"- {item}")
        return 1
    print("KB lint passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
