# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Checks the HTML escaping helper shared by the dashboard and Genome Bins scripts."""

import json
import os
import re
import shutil
import subprocess

import pytest

JS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "js")


def _escape_with(script_name: str, value: str) -> str:
    with open(os.path.join(JS_DIR, script_name), encoding="utf-8") as handle:
        source = handle.read()
    match = re.search(r"function _escapeHtml\(s\) \{.*?\n    \}", source, re.S)
    assert match, script_name
    program = match.group(0) + f"\nprocess.stdout.write(_escapeHtml({json.dumps(value)}));"
    return subprocess.run(["node", "-e", program], capture_output=True, text=True, check=True).stdout


@pytest.mark.skipif(shutil.which("node") is None, reason="node not available")
@pytest.mark.parametrize("script_name", ["dashboard.js", "bins.js"])
def test_escape_html_handles_attribute_quotes(script_name):
    assert _escape_with(script_name, "a\"b'<c>&") == "a&quot;b&#39;&lt;c&gt;&amp;"
    assert _escape_with(script_name, "") == ""


@pytest.mark.skipif(shutil.which("node") is None, reason="node not available")
def test_bin_names_break_after_separators_and_stay_escaped():
    with open(os.path.join(JS_DIR, "bins.js"), encoding="utf-8") as handle:
        source = handle.read()
    helpers = [re.search(r"function %s\(s\) \{.*?\n    \}" % name, source, re.S)
               for name in ("_escapeHtml", "_breakableName")]
    assert all(helpers)
    program = "\n".join(match.group(0) for match in helpers) + (
        "\nprocess.stdout.write(JSON.stringify([_breakableName('output/bin_1.a-b'),"
        " _breakableName('<b>&x')]));")
    result = subprocess.run(["node", "-e", program], capture_output=True, text=True, check=True).stdout
    assert json.loads(result) == ["output/<wbr>bin_<wbr>1.<wbr>a-<wbr>b", "&lt;b&gt;&amp;x"]
