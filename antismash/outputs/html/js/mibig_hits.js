/* Standalone MIBiG-hits page renderer.
 *
 * Loaded by mibig_hits.html, which a gene's "MIBiG Hits" link opens in a new tab as
 *     mibig_hits.html?anchor=<anchor>&gene=<gene>
 * This is a real file:// page (one shared file, not one per gene), so it loads jQuery +
 * tablesorter + pako via relative paths and renders exactly like the upstream per-gene page.
 *
 * The per-region hit data lives in knownclusterblast/<anchor>_hits.js as
 *     mibigHits[anchor] = "<base64 of gzipped JSON>"
 * decoded to {gene: [[8 columns], ...]}. It is injected on demand (works under file://),
 * decompressed with pako, and rendered as the sortable table.
 */
"use strict";

var mibigHits = mibigHits || {};   // anchor -> base64(gzip(JSON)), filled by the injected data file

var _MIBIG_COLUMNS = ["MIBiG Protein", "Description", "MIBiG Cluster", "MiBiG Product",
                      "% ID", "% Coverage", "BLAST Score", "E-value"];

function _mibigEsc(value) {
    var div = document.createElement("div");
    div.textContent = String(value === undefined || value === null ? "" : value);
    return div.innerHTML;
}

function _mibigParams() {
    var params = {};
    (window.location.search || "").replace(/^\?/, "").split("&").forEach(function (pair) {
        if (!pair) { return; }
        var parts = pair.split("=");
        params[decodeURIComponent(parts[0])] = decodeURIComponent((parts[1] || "").replace(/\+/g, " "));
    });
    return params;
}

function _mibigDecode(anchor) {
    var raw = mibigHits[anchor];
    if (!raw) { return null; }
    var bytes = Uint8Array.from(atob(raw), function (c) { return c.charCodeAt(0); });
    return JSON.parse(pako.inflate(bytes, { to: "string" }));
}

function _mibigRenderTable(rows) {
    var html = '<table id="myTable" class="tablesorter" border="1"><thead><tr>';
    _MIBIG_COLUMNS.forEach(function (col) { html += "<th>" + _mibigEsc(col) + "</th>"; });
    html += "</tr></thead><tbody>";
    rows.forEach(function (row) {
        html += "<tr>" +
            "<td>" + _mibigEsc(row[0]) + "</td>" +
            "<td>" + _mibigEsc(row[1]) + "</td>" +
            '<td><a target="_blank" rel="noopener" href="https://mibig.secondarymetabolites.org/go/' +
                encodeURIComponent(row[2]) + '">' + _mibigEsc(row[2]) + "</a></td>" +
            "<td>" + _mibigEsc(row[3]) + "</td>" +
            "<td>" + _mibigEsc(row[4]) + "</td>" +
            "<td>" + _mibigEsc(row[5]) + "</td>" +
            "<td>" + _mibigEsc(row[6]) + "</td>" +
            "<td>" + _mibigEsc(row[7]) + "</td>" +
            "</tr>";
    });
    html += "</tbody></table>";
    return html;
}

function _mibigDraw() {
    var params = _mibigParams();
    var anchor = params.anchor;
    var gene = params.gene;
    var container = document.getElementById("mibig-content") || document.body;
    var rows = [];
    try {
        var data = _mibigDecode(anchor) || {};
        rows = data[gene] || [];
    } catch (err) {
        container.innerHTML = "<p>Could not read MIBiG hit data.</p>";
        return;
    }
    document.title = gene + " MIBiG hits";
    container.innerHTML = "<h3>" + _mibigEsc(gene) + " &mdash; MIBiG hits (" + rows.length + ")</h3>" +
        (rows.length ? _mibigRenderTable(rows) : "<p>No MIBiG hits for this gene.</p>");
    if (rows.length && window.jQuery && jQuery.fn && jQuery.fn.tablesorter) {
        jQuery("#myTable").tablesorter();
    }
}

(function () {
    function start() {
        var anchor = _mibigParams().anchor;
        if (!anchor) {
            (document.getElementById("mibig-content") || document.body).innerHTML =
                "<p>No region specified.</p>";
            return;
        }
        if (mibigHits[anchor]) { _mibigDraw(); return; }
        var script = document.createElement("script");
        script.src = "knownclusterblast/" + anchor + "_hits.js";
        script.onload = _mibigDraw;
        script.onerror = function () {
            (document.getElementById("mibig-content") || document.body).innerHTML =
                "<p>MIBiG hit data could not be loaded.</p>";
        };
        document.head.appendChild(script);
    }
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", start);
    } else {
        start();
    }
})();
