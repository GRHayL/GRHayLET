#!/usr/bin/env python3
"""Deterministic IllinoisGRMHD knowledge-base checker."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
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
    r"(?<![\w/])IllinoisGRMHD/[A-Za-z0-9_./*?+\[\]-]+"
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
    files.extend(sorted(WIKI.rglob("*.md")))
    return sorted({path.resolve() for path in files})


def wiki_pages() -> list[Path]:
    """Return live wiki Markdown pages."""
    return sorted(path.resolve() for path in WIKI.rglob("*.md"))


def mask_code(text: str) -> str:
    """Mask fenced code and HTML comments while preserving line count."""
    pattern = re.compile(
        r"^[ ]{0,3}(?P<fence>`{3,}|~{3,})[^\n]*\n.*?^[ ]{0,3}(?P=fence)[ \t]*$|<!--.*?-->",
        re.MULTILINE | re.DOTALL,
    )

    def repl(match: re.Match[str]) -> str:
        return "".join("\n" if char == "\n" else " " for char in match.group(0))

    return pattern.sub(repl, text)


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


def source_path_mentions(text: str) -> list[tuple[int, str]]:
    """Return IllinoisGRMHD path/glob mentions independent of Markdown style."""
    masked = mask_code(text)
    return [
        (masked.count("\n", 0, match.start()) + 1, match.group(0).rstrip(".,;:"))
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
    live = {path.relative_to(WIKI).as_posix() for path in wiki_pages()}
    for missing in sorted(EXPECTED_WIKI - live):
        fail(failures, WIKI / missing, "required wiki page missing")
    for extra in sorted(live - EXPECTED_WIKI):
        fail(failures, WIKI / extra, "unexpected wiki page outside exact tree")
    raw_files = {
        path.relative_to(ROOT / "raw").as_posix()
        for path in (ROOT / "raw").rglob("*")
        if path.is_file()
    }
    for extra in sorted(raw_files - {"SOURCES.md"}):
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
    allowed_global = {path.resolve() for path in SUPPORT_PAGES} | {AGENTS.resolve(), SOURCES.resolve()}
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
        for missing in sorted(child_paths - targets):
            fail(failures, router, f"router missing immediate child: {rel(missing)}")
        for target in sorted(targets):
            if target not in child_paths and target not in allowed_global:
                fail(failures, router, f"illegal router edge: {rel(target)}")


def check_leaves(failures: list[str]) -> None:
    """Check leaf header, section, parent, source, and neighbor contract."""
    for path in wiki_pages():
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
    live = set(wiki_pages())
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

    if failures:
        print("KB lint failed:")
        for item in sorted(set(failures)):
            print(f"- {item}")
        return 1
    print("KB lint passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
