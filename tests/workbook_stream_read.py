"""Read-only OOXML/streaming acceptance. Never repairs or writes an XLSX."""
import argparse
import hashlib
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile
import openpyxl
from openpyxl.utils.cell import coordinate_to_tuple, get_column_letter


def actual_range(archive, member):
    lo_r = lo_c = None
    hi_r = hi_c = 0
    with archive.open(member) as handle:
        for _, node in ET.iterparse(handle, events=("end",)):
            if node.tag.rsplit("}", 1)[-1] == "c":
                populated = any(child.tag.rsplit("}", 1)[-1] in ("v", "is", "f") for child in node)
                if populated:
                    r, c = coordinate_to_tuple(node.attrib["r"])
                    lo_r, lo_c = min(lo_r or r, r), min(lo_c or c, c)
                    hi_r, hi_c = max(hi_r, r), max(hi_c, c)
                node.clear()
            elif node.tag.rsplit("}", 1)[-1] == "row":
                node.clear()
    if not hi_r:
        return "A1"
    return f"{get_column_letter(lo_c)}{lo_r}:{get_column_letter(hi_c)}{hi_r}"


def matches(column, i, cell):
    expected = column["values"][i]
    actual = cell.value
    if column["missing"][i]:
        return actual == "NA" and cell.data_type in ("s", "inlineStr")
    kind = column["kind"]
    if kind == "numeric":
        if not isinstance(actual, (int, float)) or isinstance(actual, bool):
            return False
        expected = float(expected)
        if not math.isfinite(expected) or not math.isfinite(actual):
            return expected == actual
        equal = abs(expected - actual) <= 1e-12 * max(1.0, abs(expected))
        if re.search(r"(^p$|pval|qval|formal_q|q_bh|padj|p[.]adj|BH|p_value|p[.]value|Pr\()", column["name"]):
            equal = equal and ((expected < .05) == (actual < .05)) and ((expected == 0) == (actual == 0))
        return equal
    if kind == "logical":
        return isinstance(actual, bool) and actual == expected
    if expected == "" and actual is None:
        return True
    return isinstance(actual, str) and actual.replace("\r\n", "\n") == str(expected).replace("\r\n", "\n")


def validate(workbook, specification):
    # read_only uses the declared dimensions by default. Do not override them.
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=False)
    results = []
    try:
        expected_sheets = specification["sheets"]
        if book.sheetnames != [s["name"] for s in expected_sheets]:
            raise ValueError("Sheet names/order differ from the expected written tables")
        with ZipFile(workbook) as archive:
            for spec in expected_sheets:
                ws = book[spec["name"]]
                declared = ws.calculate_dimension()
                # _worksheet_path is supplied by openpyxl's workbook relationship parser.
                occupied = actual_range(archive, ws._worksheet_path)
                row = dict(sheet=ws.title, declared_range=declared, actual_range=occupied,
                           expected_range=spec["expected_range"], expected_rows=spec["rows"],
                           expected_columns=len(spec["columns"]), read_rows=0, read_columns=0,
                           cell_type_and_value_pass=True, pass_status=False, errors=[])
                if declared != spec["expected_range"] or occupied != spec["expected_range"]:
                    row["errors"].append("Declared/actual/expected worksheet ranges differ")
                for r, cells in enumerate(ws.iter_rows(), 1):
                    row["read_rows"] = r
                    row["read_columns"] = max(row["read_columns"], len(cells))
                    if r > spec["rows"] or len(cells) != len(spec["columns"]):
                        row["errors"].append(f"Unexpected row/column extent at row {r}")
                        break
                    if r == 1:
                        if [c.value for c in cells] != [c["name"] for c in spec["columns"]]:
                            row["errors"].append("Column names/order differ")
                        continue
                    for c, cell in zip(spec["columns"], cells):
                        if not matches(c, r - 2, cell):
                            row["cell_type_and_value_pass"] = False
                            if len(row["errors"]) < 30:
                                row["errors"].append(f"Type/value/missing mismatch at {cell.coordinate}: {c['name']}")
                if row["read_rows"] != spec["rows"] or row["read_columns"] != len(spec["columns"]):
                    row["errors"].append("Streaming read is incomplete")
                row["pass_status"] = not row["errors"]
                results.append(row)
    finally:
        book.close()
    return dict(reader="openpyxl", version=openpyxl.__version__, read_only=True,
                data_only=False, sheets=results, pass_status=all(r["pass_status"] for r in results))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--workbook")
    parser.add_argument("--spec")
    parser.add_argument("--report")
    args = parser.parse_args()
    if args.probe:
        print(json.dumps(dict(python=sys.version.split()[0], openpyxl=openpyxl.__version__, read_only=True)))
        return 0
    if not all((args.workbook, args.spec, args.report)):
        parser.error("--workbook, --spec and --report are required")
    try:
        result = validate(args.workbook, json.loads(Path(args.spec).read_text(encoding="utf-8")))
    except Exception as exc:
        result = dict(pass_status=False, error=f"{type(exc).__name__}: {exc}", sheets=[])
    with open(args.workbook, "rb") as handle:
        result["workbook_sha256"] = hashlib.file_digest(handle, "sha256").hexdigest()
    Path(args.report).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print("PASS streaming workbook read" if result["pass_status"] else "FAIL streaming workbook read; see " + args.report)
    return 0 if result["pass_status"] else 1


if __name__ == "__main__":
    sys.exit(main())
