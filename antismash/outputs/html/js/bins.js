/* Genome Bins page logic for MetaSMASH — sortable, filterable bin table */

$(document).ready(function () {
    if (typeof binsData === "undefined" || !binsData) {
        return;
    }

    // The selection is a one-shot handoff to regions.html; returning here
    // consumes it so the next visit to regions.html shows all regions.
    sessionStorage.removeItem("antismash_bgc_selection");
    sessionStorage.removeItem("metasmash_bgc_selection_label");

    var UNBINNED_SELECTION_LABEL = "unbinned contigs";

    // Server order is BGCs descending, then bin name
    var _bins = (binsData.bins || []).slice();
    var _unbinned = binsData.unbinned;
    var _sortColumn = null;
    var _sortAsc = true;
    var _numericKeys = {contigs: true, bgcs: true, complete: true, on_edge: true, known_hits: true};

    function _escapeHtml(s) {
        if (!s) return "";
        // safe for both element content and quoted attribute values
        return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
            .replace(/"/g, "&quot;").replace(/'/g, "&#39;");
    }

    // Escaped name with line-break hints after separators, so wrapping
    // keeps "cleanbin_000067" readable rather than splitting mid-word
    function _breakableName(s) {
        return _escapeHtml(s).replace(/([\/_.\-])/g, "$1<wbr>");
    }

    function _formatNumber(n) {
        return n.toLocaleString();
    }

    function _categoryLabel(pair) {
        return pair[1] > 1 ? pair[0] + " ×" + pair[1] : pair[0];
    }

    // Category colours come from the region button classes, as on the dashboard chart
    var _colorCache = {};
    function _categoryColor(category) {
        if (_colorCache[category]) return _colorCache[category];
        var el = document.createElement("div");
        el.className = "regbutton " + category;
        el.style.position = "absolute";
        el.style.visibility = "hidden";
        document.body.appendChild(el);
        var bg = window.getComputedStyle(el).backgroundColor;
        document.body.removeChild(el);
        if (!bg || bg === "transparent" || bg === "rgba(0, 0, 0, 0)" || bg === "rgb(255, 255, 255)") {
            bg = "#888";
        }
        _colorCache[category] = bg;
        return bg;
    }

    // White text on dark category colours, black on light ones
    function _textColor(background) {
        var rgb = background.match(/\d+/g);
        if (!rgb || rgb.length < 3) return "#000";
        var luminance = (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) / 255;
        return luminance < 0.5 ? "#fff" : "#000";
    }

    function renderCategories(categories) {
        if (!categories.length) return "-";
        return categories.map(function (pair) {
            var background = _categoryColor(pair[0]);
            return '<span class="bgc-pill bgc-category-pill" style="background-color: '
                + _escapeHtml(background) + '; color: ' + _textColor(background) + '">'
                + _escapeHtml(_categoryLabel(pair)) + '</span>';
        }).join("");
    }

    function renderRow(bin, index, isUnbinned) {
        var cls = isUnbinned ? "bins-unbinned-row" : (index % 2 === 0 ? "odd" : "even");
        var html = ['<tr class="' + cls + '">'];
        html.push('<td class="bins-name-cell">' + _breakableName(bin.bin) + '</td>');
        html.push('<td class="digits">' + (bin.contigs === null ? "–" : _formatNumber(bin.contigs)) + '</td>');
        html.push('<td class="digits">' + _formatNumber(bin.bgcs) + '</td>');
        html.push('<td class="digits">' + _formatNumber(bin.complete) + '</td>');
        html.push('<td class="digits">' + _formatNumber(bin.on_edge) + '</td>');
        html.push('<td class="digits">' + _formatNumber(bin.known_hits) + '</td>');
        html.push('<td class="bins-categories-cell">' + renderCategories(bin.categories) + '</td>');
        if (bin.bgcs > 0) {
            var target = isUnbinned ? "unbinned" : String(index);
            html.push('<td><a href="regions.html" class="bins-view-link" data-target="' + target + '">View BGCs &rarr;</a></td>');
        } else {
            html.push('<td></td>');
        }
        html.push('</tr>');
        return html.join("");
    }

    function render() {
        var query = ($("#bins-filter").val() || "").toLowerCase();
        var showEmpty = $("#bins-show-empty").prop("checked");
        var html = [];
        var shown = 0;
        for (var i = 0; i < _bins.length; i++) {
            var bin = _bins[i];
            if (!showEmpty && bin.bgcs === 0) continue;
            if (query && bin.bin.toLowerCase().indexOf(query) === -1) continue;
            html.push(renderRow(bin, i, false));
            shown++;
        }
        // The unbinned row always stays last, whatever the sort or filter
        if (_unbinned) {
            html.push(renderRow(_unbinned, shown, true));
        }
        document.getElementById("bins-table-body").innerHTML = html.join("");
        $("#bins-table-title").text("Genome Bins (" + shown + " of " + _bins.length + ")");
    }

    function sortBins(col, asc) {
        _sortColumn = col;
        _sortAsc = asc;
        var isNumeric = _numericKeys[col];
        _bins.sort(function (a, b) {
            var va = a[col];
            var vb = b[col];
            if (isNumeric) {
                return asc ? va - vb : vb - va;
            }
            va = String(va).toLowerCase();
            vb = String(vb).toLowerCase();
            if (va < vb) return asc ? -1 : 1;
            if (va > vb) return asc ? 1 : -1;
            return 0;
        });
        render();
    }

    $("#bins-table thead th").on("click", function () {
        if ($(this).attr("data-no-sort") === "true") return;
        var key = $(this).attr("data-sort-key");
        if (!key) return;
        var asc = (_sortColumn === key) ? !_sortAsc : true;
        $("#bins-table thead th").removeClass("sort-asc sort-desc");
        $(this).addClass(asc ? "sort-asc" : "sort-desc");
        sortBins(key, asc);
    });

    var _filterTimer = null;
    $("#bins-filter").on("keyup", function () {
        clearTimeout(_filterTimer);
        _filterTimer = setTimeout(render, 150);
    });

    $("#bins-show-empty").on("change", render);

    // "View BGCs" — hand the bin's regions to regions.html, as the dashboard's
    // "View detailed results for selection" does
    $("#bins-table-body").on("click", ".bins-view-link", function (e) {
        e.preventDefault();
        var target = $(this).attr("data-target");
        var bin = target === "unbinned" ? _unbinned : _bins[parseInt(target, 10)];
        if (!bin || !bin.anchors.length) return;
        try {
            sessionStorage.setItem("antismash_bgc_selection", JSON.stringify(bin.anchors));
            sessionStorage.setItem("metasmash_bgc_selection_label",
                                   target === "unbinned" ? UNBINNED_SELECTION_LABEL : bin.bin);
        } catch (err) {
            console.warn("Could not store the bin selection", err);
        }
        window.location.href = "regions.html";
    });

    render();
});
