#!/usr/bin/env python3
"""Validate and repair the OOXML package structure used by Promotheus outputs."""

from __future__ import annotations

import os
import posixpath
import re
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree as ET


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CONTENT_TYPES_NS = "http://schemas.openxmlformats.org/package/2006/content-types"

DRAWING_RELATIONSHIPS = {
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships/drawing",
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships/vmlDrawing",
}


def _source_part_for_relationships(relationships_part: str) -> str:
    path = PurePosixPath(relationships_part)
    if relationships_part == "_rels/.rels":
        return ""
    if path.parent.name != "_rels" or not path.name.endswith(".rels"):
        raise ValueError(f"invalid relationships part path: {relationships_part}")
    return str(path.parent.parent / path.name.removesuffix(".rels"))


def _relationship_target(relationships_part: str, target: str) -> str:
    if target.startswith("/"):
        return target.lstrip("/")
    source_part = _source_part_for_relationships(relationships_part)
    source_dir = posixpath.dirname(source_part)
    return posixpath.normpath(posixpath.join(source_dir, target)).lstrip("/")


def _cell_coordinates(reference: str) -> tuple[int, int]:
    match = re.fullmatch(r"\$?([A-Za-z]+)\$?([1-9][0-9]*)", reference)
    if match is None:
        raise ValueError(f"invalid cell reference: {reference}")
    column = 0
    for character in match.group(1).upper():
        column = column * 26 + ord(character) - ord("A") + 1
    return int(match.group(2)), column


def _column_name(column: int) -> str:
    name = ""
    while column:
        column, remainder = divmod(column - 1, 26)
        name = chr(ord("A") + remainder) + name
    return name


def _cell_bounds(worksheet_root: ET.Element) -> tuple[int, int, int, int] | None:
    coordinates = []
    for cell in worksheet_root.findall(f".//{{{MAIN_NS}}}c"):
        reference = cell.attrib.get("r")
        if reference:
            coordinates.append(_cell_coordinates(reference))
    if not coordinates:
        return None
    rows = [row for row, _ in coordinates]
    columns = [column for _, column in coordinates]
    return min(rows), min(columns), max(rows), max(columns)


def _bounds_reference(bounds: tuple[int, int, int, int]) -> str:
    min_row, min_column, max_row, max_column = bounds
    start = f"{_column_name(min_column)}{min_row}"
    end = f"{_column_name(max_column)}{max_row}"
    return start if start == end else f"{start}:{end}"


def _dimension_covers(dimension: str, bounds: tuple[int, int, int, int]) -> bool:
    start, _, end = dimension.replace("$", "").partition(":")
    if not end:
        end = start
    dim_min_row, dim_min_column = _cell_coordinates(start)
    dim_max_row, dim_max_column = _cell_coordinates(end)
    min_row, min_column, max_row, max_column = bounds
    return (
        dim_min_row <= min_row
        and dim_min_column <= min_column
        and dim_max_row >= max_row
        and dim_max_column >= max_column
    )


def validate_xlsx_package(path: Path) -> list[str]:
    """Return structural OOXML errors without changing *path*."""
    errors: list[str] = []
    try:
        with zipfile.ZipFile(path) as archive:
            bad_member = archive.testzip()
            if bad_member:
                errors.append(f"corrupt ZIP member: {bad_member}")
            members = set(archive.namelist())
            required = {"[Content_Types].xml", "_rels/.rels", "xl/workbook.xml", "xl/_rels/workbook.xml.rels"}
            for member in sorted(required - members):
                errors.append(f"missing required XLSX part: {member}")

            for relationships_part in sorted(member for member in members if member.endswith(".rels")):
                try:
                    root = ET.fromstring(archive.read(relationships_part))
                except ET.ParseError as exc:
                    errors.append(f"invalid relationships XML {relationships_part}: {exc}")
                    continue
                for relationship in root.findall(f"{{{PKG_REL_NS}}}Relationship"):
                    if relationship.attrib.get("TargetMode") == "External":
                        continue
                    target = relationship.attrib.get("Target", "")
                    if not target:
                        errors.append(f"empty relationship target in {relationships_part}")
                        continue
                    resolved = _relationship_target(relationships_part, target)
                    if resolved not in members:
                        errors.append(f"dangling relationship in {relationships_part}: {resolved}")

            if "[Content_Types].xml" in members:
                try:
                    content_types = ET.fromstring(archive.read("[Content_Types].xml"))
                    for override in content_types.findall(f"{{{CONTENT_TYPES_NS}}}Override"):
                        part = override.attrib.get("PartName", "").lstrip("/")
                        if part and part not in members:
                            errors.append(f"content type references missing part: {part}")
                except ET.ParseError as exc:
                    errors.append(f"invalid [Content_Types].xml: {exc}")

            for worksheet_part in sorted(
                member
                for member in members
                if re.fullmatch(r"xl/worksheets/[^/]+\.xml", member)
            ):
                try:
                    worksheet = ET.fromstring(archive.read(worksheet_part))
                    bounds = _cell_bounds(worksheet)
                    dimension = worksheet.find(f"{{{MAIN_NS}}}dimension")
                    if bounds is not None and dimension is not None:
                        reference = dimension.attrib.get("ref", "")
                        try:
                            covers = bool(reference) and _dimension_covers(reference, bounds)
                        except ValueError:
                            covers = False
                        if not covers:
                            errors.append(
                                f"worksheet dimension {reference or '<empty>'} in {worksheet_part} "
                                f"does not cover populated range {_bounds_reference(bounds)}"
                            )
                except (ET.ParseError, ValueError) as exc:
                    errors.append(f"invalid worksheet XML {worksheet_part}: {exc}")
    except (OSError, zipfile.BadZipFile) as exc:
        errors.append(f"invalid XLSX ZIP package: {exc}")
    return errors


def repair_xlsx_package(path: Path) -> list[str]:
    """Repair known openxlsx 4.2.5.2 packaging defects atomically."""
    with zipfile.ZipFile(path) as archive:
        members = set(archive.namelist())
        entries = [(info, archive.read(info.filename)) for info in archive.infolist()]

    repaired: list[str] = []
    rewritten: dict[str, bytes] = {}
    for info, data in entries:
        name = info.filename
        if name.endswith(".rels"):
            root = ET.fromstring(data)
            text = data.decode("utf-8")
            for relationship in list(root.findall(f"{{{PKG_REL_NS}}}Relationship")):
                if relationship.attrib.get("TargetMode") == "External":
                    continue
                target = relationship.attrib.get("Target", "")
                resolved = _relationship_target(name, target) if target else ""
                if resolved not in members and relationship.attrib.get("Type") in DRAWING_RELATIONSHIPS:
                    relationship_id = re.escape(relationship.attrib.get("Id", ""))
                    pattern = rf'<(?:[A-Za-z_][\w.-]*:)?Relationship\b(?=[^>]*\bId=["\']{relationship_id}["\'])[^>]*/>'
                    text, count = re.subn(pattern, "", text, count=1)
                    if count != 1:
                        raise ValueError(f"cannot remove relationship {relationship.attrib.get('Id')} from {name}")
                    repaired.append(f"removed dangling drawing relationship {resolved}")
            if text != data.decode("utf-8"):
                rewritten[name] = text.encode("utf-8")
        elif name == "[Content_Types].xml":
            root = ET.fromstring(data)
            text = data.decode("utf-8")
            for override in list(root.findall(f"{{{CONTENT_TYPES_NS}}}Override")):
                part = override.attrib.get("PartName", "").lstrip("/")
                if part and part not in members and part.startswith("xl/drawings/"):
                    part_name = re.escape(override.attrib.get("PartName", ""))
                    pattern = rf'<(?:[A-Za-z_][\w.-]*:)?Override\b(?=[^>]*\bPartName=["\']{part_name}["\'])[^>]*/>'
                    text, count = re.subn(pattern, "", text, count=1)
                    if count != 1:
                        raise ValueError(f"cannot remove content type override for {part}")
                    repaired.append(f"removed content type for missing drawing {part}")
            if text != data.decode("utf-8"):
                rewritten[name] = text.encode("utf-8")
        elif re.fullmatch(r"xl/worksheets/[^/]+\.xml", name):
            root = ET.fromstring(data)
            bounds = _cell_bounds(root)
            dimension = root.find(f"{{{MAIN_NS}}}dimension")
            if bounds is not None and dimension is not None:
                reference = dimension.attrib.get("ref", "")
                try:
                    covers = bool(reference) and _dimension_covers(reference, bounds)
                except ValueError:
                    covers = False
                if not covers:
                    text = data.decode("utf-8")
                    dimension_tag = re.search(r'<(?:[A-Za-z_][\w.-]*:)?dimension\b[^>]*/?>', text)
                    if dimension_tag is None:
                        raise ValueError(f"cannot locate worksheet dimension in {name}")
                    replacement = re.sub(
                        r'\bref=(["\']).*?\1',
                        f'ref="{_bounds_reference(bounds)}"',
                        dimension_tag.group(0),
                        count=1,
                    )
                    text = text[: dimension_tag.start()] + replacement + text[dimension_tag.end() :]
                    rewritten[name] = text.encode("utf-8")
                    repaired.append(f"corrected worksheet dimension in {name}")

    if not repaired:
        return repaired

    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, delete=False) as handle:
            temporary_name = handle.name
        with zipfile.ZipFile(temporary_name, "w") as output:
            for info, data in entries:
                output.writestr(info, rewritten.get(info.filename, data))
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)
    return repaired
