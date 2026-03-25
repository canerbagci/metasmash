/* Dashboard page logic for MetaSMASH — client-side paginated table */

$(document).ready(function () {
    if (typeof dashboardData === "undefined" || !dashboardData) {
        return;
    }

    // =========================================================================
    //  Table state
    // =========================================================================
    var _allRows = dashboardData.bgc_table || [];
    var _filteredRows = _allRows;
    var _currentPage = 0;
    var _pageSize = 100;
    var _sortColumn = null;
    var _sortAsc = true;
    var _selectedAnchors = {};  // anchor_id -> true
    var _hasTaxonomy = dashboardData.has_taxonomy || false;
    var _meta = dashboardData.filter_metadata || {};

    // Pre-compute a sortable region label "record_index.region_number"
    for (var i = 0; i < _allRows.length; i++) {
        _allRows[i].region_label = _allRows[i].record_index + "." + _allRows[i].region_number;
    }

    // =========================================================================
    //  Centralized filter state
    // =========================================================================
    var _filterState = {
        text: "",
        chartLabel: null,
        chartViewMode: null,
        edgeMode: "all",         // "all" | "complete" | "on_edge"
        similarityMin: _meta.similarity_min || 0,
        similarityMax: _meta.similarity_max || 100,
        selectedTypes: [],       // empty = all
        selectedCategories: []   // empty = all
    };

    // Store defaults so we can detect active filters
    var _defaults = {
        similarityMin: _filterState.similarityMin,
        similarityMax: _filterState.similarityMax
    };

    function applyAllFilters() {
        var text = _filterState.text ? _filterState.text.toLowerCase() : "";
        var chartLabel = _filterState.chartLabel;
        var chartKey = _filterState.chartViewMode === "type" ? "products" : "product_categories";
        var edge = _filterState.edgeMode;
        var sMin = _filterState.similarityMin;
        var sMax = _filterState.similarityMax;
        var types = _filterState.selectedTypes;
        var cats = _filterState.selectedCategories;
        var simActive = sMin !== _defaults.similarityMin || sMax !== _defaults.similarityMax;

        _filteredRows = _allRows.filter(function (row) {
            // Text search
            if (text) {
                var match = (row.record_id && row.record_id.toLowerCase().indexOf(text) > -1)
                    || (row.products && row.products.toLowerCase().indexOf(text) > -1)
                    || (row.product_categories && row.product_categories.toLowerCase().indexOf(text) > -1)
                    || (row.most_similar_known && row.most_similar_known.toLowerCase().indexOf(text) > -1)
                    || (row.region_label && row.region_label.indexOf(text) > -1)
                    || (_hasTaxonomy && row.taxonomy && row.taxonomy.toLowerCase().indexOf(text) > -1);
                if (!match) return false;
            }

            // Chart label filter
            if (chartLabel) {
                if (!row[chartKey] || row[chartKey].split(", ").indexOf(chartLabel) === -1) {
                    return false;
                }
            }

            // Edge filter
            if (edge === "complete" && row.contig_edge) return false;
            if (edge === "on_edge" && !row.contig_edge) return false;

            // Similarity range
            if (simActive) {
                if (row.similarity_percentage < 0) return false;
                if (row.similarity_percentage < sMin || row.similarity_percentage > sMax) return false;
            }

            // Type multi-select
            if (types.length > 0) {
                var rowTypes = row.products ? row.products.split(", ") : [];
                var typeMatch = false;
                for (var t = 0; t < types.length; t++) {
                    if (rowTypes.indexOf(types[t]) > -1) { typeMatch = true; break; }
                }
                if (!typeMatch) return false;
            }

            // Category multi-select
            if (cats.length > 0) {
                var rowCats = row.product_categories ? row.product_categories.split(", ") : [];
                var catMatch = false;
                for (var c = 0; c < cats.length; c++) {
                    if (rowCats.indexOf(cats[c]) > -1) { catMatch = true; break; }
                }
                if (!catMatch) return false;
            }

            return true;
        });

        if (_sortColumn) {
            sortRows(_sortColumn, _sortAsc);
        } else {
            _currentPage = 0;
            renderPage();
        }

        renderChips();
    }

    // =========================================================================
    //  Rendering
    // =========================================================================
    var _tbody = document.getElementById("bgc-table-body");

    function _formatNumber(n) {
        return n.toLocaleString();
    }

    function _escapeHtml(s) {
        if (!s) return "";
        var div = document.createElement("div");
        div.appendChild(document.createTextNode(s));
        return div.innerHTML;
    }

    function renderPage() {
        var start = _currentPage * _pageSize;
        var end = Math.min(start + _pageSize, _filteredRows.length);
        var html = [];

        for (var i = start; i < end; i++) {
            var bgc = _filteredRows[i];
            var cycle = (i - start) % 2 === 0 ? "odd" : "even";
            var checked = _selectedAnchors[bgc.anchor_id] ? " checked" : "";

            html.push('<tr class="' + cycle + '" data-anchor="' + bgc.anchor_id + '">');
            html.push('<td class="bgc-select-cell"><input type="checkbox" class="bgc-select-checkbox" data-anchor="' + bgc.anchor_id + '"' + checked + ' /></td>');
            html.push('<td>' + _escapeHtml(bgc.record_id) + '</td>');
            if (_hasTaxonomy) {
                html.push('<td class="taxonomy-cell" data-tooltip="' + _escapeHtml(bgc.taxonomy) + '">' + _escapeHtml(bgc.taxonomy) + '</td>');
            }
            html.push('<td>Region ' + bgc.record_index + '.' + bgc.region_number + '</td>');
            html.push('<td>' + _escapeHtml(bgc.products) + '</td>');
            html.push('<td>' + _escapeHtml(bgc.product_categories) + '</td>');
            html.push('<td class="digits">' + _formatNumber(bgc.start) + '</td>');
            html.push('<td class="digits">' + _formatNumber(bgc.end) + '</td>');
            html.push('<td class="digits">' + _formatNumber(bgc.length) + '</td>');
            html.push('<td class="digits">' + bgc.gene_count + '</td>');

            if (bgc.contig_edge) {
                html.push('<td><span class="edge-badge-yes">Yes</span></td>');
            } else {
                html.push('<td><span class="edge-badge-no">No</span></td>');
            }

            // Most similar known cluster
            var simCell = '-';
            if (bgc.most_similar_known) {
                var simClass = "similarity-low";
                if (bgc.similarity_percentage > 75) simClass = "similarity-high";
                else if (bgc.similarity_percentage > 50) simClass = "similarity-medium";
                var simContent;
                if (bgc.most_similar_url) {
                    simContent = '<a href="' + _escapeHtml(bgc.most_similar_url) + '" target="_blank">' + _escapeHtml(bgc.most_similar_known) + '</a>';
                } else {
                    simContent = _escapeHtml(bgc.most_similar_known);
                }
                simCell = '<span class="' + simClass + '" style="padding:0.1em 0.4em;border-radius:2px;">' + simContent + ' (' + bgc.similarity_percentage + '%)</span>';
            }
            html.push('<td>' + simCell + '</td>');
            html.push('<td class="bgc-navigate-cell"><a href="regions.html#' + bgc.anchor_id + '" class="bgc-navigate-link" title="View region details">&#8250;</a></td>');
            html.push('</tr>');
        }

        _tbody.innerHTML = html.join("");
        updatePaginationControls();
        updateSelectionUI();
    }

    function updatePaginationControls() {
        var totalPages = Math.max(1, Math.ceil(_filteredRows.length / _pageSize));
        var $info = $("#bgc-page-info");
        var start = _currentPage * _pageSize + 1;
        var end = Math.min((_currentPage + 1) * _pageSize, _filteredRows.length);
        if (_filteredRows.length === 0) {
            $info.text("No results");
        } else {
            $info.text(start + "-" + end + " of " + _filteredRows.length);
        }
        $("#bgc-page-prev").prop("disabled", _currentPage === 0);
        $("#bgc-page-next").prop("disabled", _currentPage >= totalPages - 1);

        // Update title count
        $("#bgc-table-title").text("All Biosynthetic Gene Clusters (" + _filteredRows.length + ")");
    }

    // =========================================================================
    //  Sorting
    // =========================================================================
    var _numericKeys = {start: 1, end: 1, length: 1, gene_count: 1, similarity_percentage: 1};

    function sortRows(col, asc) {
        _sortColumn = col;
        _sortAsc = asc;
        var isNumeric = _numericKeys[col];

        _filteredRows.sort(function (a, b) {
            var va = a[col];
            var vb = b[col];
            if (va == null) va = isNumeric ? -Infinity : "";
            if (vb == null) vb = isNumeric ? -Infinity : "";
            if (isNumeric) {
                return asc ? va - vb : vb - va;
            }
            if (typeof va === "boolean") {
                va = va ? 1 : 0;
                vb = vb ? 1 : 0;
                return asc ? va - vb : vb - va;
            }
            va = String(va).toLowerCase();
            vb = String(vb).toLowerCase();
            if (va < vb) return asc ? -1 : 1;
            if (va > vb) return asc ? 1 : -1;
            return 0;
        });

        _currentPage = 0;
        renderPage();
    }

    // Column header click sort
    $("#bgc-table thead th").on("click", function () {
        if ($(this).attr("data-no-sort") === "true") return;
        var key = $(this).attr("data-sort-key");
        if (!key) return;

        // Toggle direction if same column
        var asc = (_sortColumn === key) ? !_sortAsc : true;

        // Update sort indicators
        $("#bgc-table thead th").removeClass("sort-asc sort-desc");
        $(this).addClass(asc ? "sort-asc" : "sort-desc");

        sortRows(key, asc);
    });

    // =========================================================================
    //  Text filtering (now sets state and calls applyAllFilters)
    // =========================================================================
    var _filterTimer = null;

    $("#bgc-filter").on("keyup", function () {
        clearTimeout(_filterTimer);
        var query = $(this).val();
        _filterTimer = setTimeout(function () {
            _filterState.text = query;
            applyAllFilters();
        }, 150);
    });

    // =========================================================================
    //  Chart filter (by type/category label) — now composes with other filters
    // =========================================================================
    function applyChartFilter(label, viewMode) {
        _filterState.chartLabel = label;
        _filterState.chartViewMode = viewMode;
        applyAllFilters();
    }

    function clearChartFilter() {
        _filterState.chartLabel = null;
        _filterState.chartViewMode = null;
        $("#chart-filter-buttons .chart-filter-active").removeClass("chart-filter-active");
        $("#chart-filter-buttons .chart-filter-all").addClass("chart-filter-active");
        applyAllFilters();
    }

    // =========================================================================
    //  Pagination controls
    // =========================================================================
    $("#bgc-page-prev").on("click", function () {
        if (_currentPage > 0) {
            _currentPage--;
            renderPage();
        }
    });

    $("#bgc-page-next").on("click", function () {
        var totalPages = Math.ceil(_filteredRows.length / _pageSize);
        if (_currentPage < totalPages - 1) {
            _currentPage++;
            renderPage();
        }
    });

    $("#bgc-page-size").on("change", function () {
        _pageSize = parseInt($(this).val(), 10);
        _currentPage = 0;
        renderPage();
    });

    // =========================================================================
    //  Striped pattern for on-edge segments
    // =========================================================================
    function createStripedPattern(ctx, baseColor) {
        var size = 10;
        var patternCanvas = document.createElement("canvas");
        patternCanvas.width = size;
        patternCanvas.height = size;
        var pctx = patternCanvas.getContext("2d");

        pctx.fillStyle = baseColor;
        pctx.fillRect(0, 0, size, size);

        pctx.strokeStyle = "rgba(255,255,255,0.5)";
        pctx.lineWidth = 2;
        pctx.beginPath();
        pctx.moveTo(0, size);
        pctx.lineTo(size, 0);
        pctx.stroke();
        pctx.beginPath();
        pctx.moveTo(-size / 2, size / 2);
        pctx.lineTo(size / 2, -size / 2);
        pctx.stroke();
        pctx.beginPath();
        pctx.moveTo(size / 2, size + size / 2);
        pctx.lineTo(size + size / 2, size / 2);
        pctx.stroke();

        return ctx.createPattern(patternCanvas, "repeat");
    }

    // =========================================================================
    //  Read product colors from CSS
    // =========================================================================
    var colorCache = {};

    function getColor(label, category) {
        var key = label;
        if (colorCache[key]) {
            return colorCache[key];
        }
        var el = document.createElement("div");
        var classes = "regbutton " + label;
        if (category && category !== label) {
            classes += " " + category;
        }
        el.className = classes;
        el.style.position = "absolute";
        el.style.visibility = "hidden";
        document.body.appendChild(el);
        var bg = window.getComputedStyle(el).backgroundColor;
        document.body.removeChild(el);
        if (!bg || bg === "transparent" || bg === "rgba(0, 0, 0, 0)" || bg === "rgb(255, 255, 255)") {
            bg = "#888";
        }
        colorCache[key] = bg;
        return bg;
    }

    // =========================================================================
    //  Chart
    // =========================================================================
    var chartCanvas = document.getElementById("bgc-chart");
    var chart = null;

    function buildChart(viewMode) {
        var countsData = viewMode === "type"
            ? dashboardData.bgc_type_counts
            : dashboardData.bgc_category_counts;

        var labels = Object.keys(countsData);
        labels.sort(function (a, b) {
            var totalA = countsData[a].complete + countsData[a].on_edge;
            var totalB = countsData[b].complete + countsData[b].on_edge;
            return totalB - totalA;
        });

        var completeValues = labels.map(function (l) { return countsData[l].complete; });
        var onEdgeValues = labels.map(function (l) { return countsData[l].on_edge; });

        var ctx = chartCanvas.getContext("2d");
        var typeToCategory = dashboardData.type_to_category || {};
        var completeColors = labels.map(function (l) { return getColor(l, typeToCategory[l]); });
        var onEdgePatterns = labels.map(function (l) {
            return createStripedPattern(ctx, getColor(l, typeToCategory[l]));
        });

        if (chart) {
            chart.destroy();
        }

        var barHeight = 28;
        var minHeight = 150;
        var computedHeight = Math.max(minHeight, labels.length * barHeight + 60);
        chartCanvas.style.height = computedHeight + "px";

        chart = new Chart(ctx, {
            type: "bar",
            data: {
                labels: labels,
                datasets: [
                    {
                        label: "Complete",
                        data: completeValues,
                        backgroundColor: completeColors,
                        borderWidth: 0
                    },
                    {
                        label: "On contig edge",
                        data: onEdgeValues,
                        backgroundColor: onEdgePatterns,
                        borderColor: completeColors,
                        borderWidth: 1
                    }
                ]
            },
            options: {
                indexAxis: "y",
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    x: {
                        stacked: true,
                        beginAtZero: true,
                        title: { display: true, text: "Number of BGCs" },
                        ticks: {
                            stepSize: 1,
                            precision: 0
                        }
                    },
                    y: {
                        stacked: true
                    }
                },
                plugins: {
                    legend: {
                        position: "top"
                    },
                    tooltip: {
                        callbacks: {
                            label: function (context) {
                                return context.dataset.label + ": " + context.parsed.x;
                            }
                        }
                    }
                }
            }
        });

        buildFilterButtons(labels, countsData, viewMode);
    }

    // =========================================================================
    //  Chart filter buttons
    // =========================================================================
    function buildFilterButtons(labels, countsData, viewMode) {
        var $container = $("#chart-filter-buttons");
        $container.empty();

        var typeToCategory = dashboardData.type_to_category || {};

        var $all = $('<div class="regbutton chart-filter-all chart-filter-active"><a>All</a></div>');
        $all.on("click", function () {
            clearChartFilter();
        });
        $container.append($all);

        labels.forEach(function (label) {
            var total = countsData[label].complete + countsData[label].on_edge;
            var classes = "regbutton " + label;
            if (viewMode === "type" && typeToCategory[label]) {
                classes += " " + typeToCategory[label];
            }
            var $btn = $('<div class="' + classes + '"><a>' + label + ' (' + total + ')</a></div>');
            $btn.on("click", function () {
                if (_filterState.chartLabel === label) {
                    clearChartFilter();
                } else {
                    $container.find(".chart-filter-active").removeClass("chart-filter-active");
                    $(this).addClass("chart-filter-active");
                    applyChartFilter(label, viewMode);
                }
            });
            $container.append($btn);
        });
    }

    if (chartCanvas && dashboardData.summary.total_bgcs > 0) {
        buildChart("category");

        $("#chart-view-toggle").on("change", function () {
            _filterState.chartLabel = null;
            _filterState.chartViewMode = null;
            buildChart($(this).val());
            applyAllFilters();
        });
    }

    // =========================================================================
    //  Advanced filter panel — collapsible toggle
    // =========================================================================
    var $filterPanel = $("#adv-filter-panel");
    var $filterBody = $("#adv-filter-body");
    var $filterArrow = $("#adv-filter-arrow");

    $("#adv-filter-toggle").on("click", function (e) {
        if ($(e.target).closest("#adv-clear-all").length) return;
        $filterBody.slideToggle(200);
        $filterArrow.toggleClass("adv-arrow-open");
    });

    // =========================================================================
    //  Edge toggle
    // =========================================================================
    $("#adv-edge-toggle .adv-toggle-btn").on("click", function () {
        var mode = $(this).data("edge");
        $("#adv-edge-toggle .adv-toggle-btn").removeClass("adv-toggle-active");
        $(this).addClass("adv-toggle-active");
        _filterState.edgeMode = mode;
        applyAllFilters();
    });

    // =========================================================================
    //  Multi-select dropdowns
    // =========================================================================
    function initMultiSelect(containerSel, items, stateKey, allLabel) {
        var $container = $(containerSel);
        var $display = $container.find(".adv-ms-display");
        var $dropdown = $container.find(".adv-ms-dropdown");
        var $search = $container.find(".adv-ms-search");
        var $options = $container.find(".adv-ms-options");

        function buildOptions(filter) {
            var q = (filter || "").toLowerCase();
            var html = "";
            for (var i = 0; i < items.length; i++) {
                if (q && items[i].toLowerCase().indexOf(q) === -1) continue;
                var checked = _filterState[stateKey].indexOf(items[i]) > -1 ? " checked" : "";
                html += '<label class="adv-ms-option"><input type="checkbox" value="' + _escapeHtml(items[i]) + '"' + checked + ' /> ' + _escapeHtml(items[i]) + '</label>';
            }
            $options.html(html);
        }

        function updateDisplay() {
            var sel = _filterState[stateKey];
            if (sel.length === 0) {
                $display.text(allLabel);
            } else if (sel.length <= 2) {
                $display.text(sel.join(", "));
            } else {
                $display.text(sel.length + " selected");
            }
        }

        // Toggle dropdown
        $display.on("click", function (e) {
            e.stopPropagation();
            var isOpen = $dropdown.is(":visible");
            // Close all other dropdowns first
            $(".adv-ms-dropdown").hide();
            if (!isOpen) {
                buildOptions("");
                $dropdown.show();
                $search.val("").focus();
            }
        });

        // Search
        $search.on("input", function () {
            buildOptions($(this).val());
        });

        // Checkbox change
        $options.on("change", "input[type=checkbox]", function () {
            var val = $(this).val();
            var arr = _filterState[stateKey];
            var idx = arr.indexOf(val);
            if ($(this).prop("checked")) {
                if (idx === -1) arr.push(val);
            } else {
                if (idx > -1) arr.splice(idx, 1);
            }
            updateDisplay();
            applyAllFilters();
        });

        // Prevent dropdown from closing when clicking inside it
        $dropdown.on("click", function (e) {
            e.stopPropagation();
        });

        // Close on outside click
        $(document).on("click", function () {
            $dropdown.hide();
        });

        buildOptions("");
        updateDisplay();

        // Return a reset function
        return function () {
            _filterState[stateKey] = [];
            buildOptions("");
            updateDisplay();
        };
    }

    var resetTypeMS = initMultiSelect("#adv-ms-type", _meta.products_list || [], "selectedTypes", "All types");
    var resetCatMS = initMultiSelect("#adv-ms-category", _meta.categories_list || [], "selectedCategories", "All categories");

    // =========================================================================
    //  Chips + Clear All
    // =========================================================================
    function renderChips() {
        var $bar = $("#adv-chips-bar");
        var chips = [];

        if (_filterState.text) {
            chips.push({label: 'Search: "' + _filterState.text + '"', key: "text"});
        }
        if (_filterState.chartLabel) {
            chips.push({label: "Chart: " + _filterState.chartLabel, key: "chart"});
        }
        if (_filterState.edgeMode !== "all") {
            var edgeLabel = _filterState.edgeMode === "complete" ? "Complete only" : "On edge only";
            chips.push({label: edgeLabel, key: "edge"});
        }
        if (_filterState.similarityMin !== _defaults.similarityMin || _filterState.similarityMax !== _defaults.similarityMax) {
            chips.push({label: "Similarity: " + _filterState.similarityMin + "% – " + _filterState.similarityMax + "%", key: "similarity"});
        }
        if (_filterState.selectedTypes.length > 0) {
            chips.push({label: "Types: " + _filterState.selectedTypes.join(", "), key: "types"});
        }
        if (_filterState.selectedCategories.length > 0) {
            chips.push({label: "Categories: " + _filterState.selectedCategories.join(", "), key: "categories"});
        }
        var html = "";
        for (var i = 0; i < chips.length; i++) {
            html += '<span class="adv-chip" data-chip-key="' + chips[i].key + '">'
                + _escapeHtml(chips[i].label)
                + ' <span class="adv-chip-remove">&times;</span></span>';
        }
        $bar.html(html);

        // Show/hide clear-all button
        $("#adv-clear-all").toggle(chips.length > 0);
    }

    // Chip removal
    $("#adv-chips-bar").on("click", ".adv-chip-remove", function () {
        var key = $(this).closest(".adv-chip").data("chip-key");
        removeFilter(key);
    });

    function removeFilter(key) {
        switch (key) {
            case "text":
                _filterState.text = "";
                $("#bgc-filter").val("");
                break;
            case "chart":
                clearChartFilter();
                return; // clearChartFilter already calls applyAllFilters
            case "edge":
                _filterState.edgeMode = "all";
                $("#adv-edge-toggle .adv-toggle-btn").removeClass("adv-toggle-active");
                $('#adv-edge-toggle [data-edge="all"]').addClass("adv-toggle-active");
                break;
            case "similarity":
                _filterState.similarityMin = _defaults.similarityMin;
                _filterState.similarityMax = _defaults.similarityMax;
                break;
            case "types":
                if (resetTypeMS) resetTypeMS();
                break;
            case "categories":
                if (resetCatMS) resetCatMS();
                break;
        }
        applyAllFilters();
    }

    // Clear all filters
    $("#adv-clear-all").on("click", function (e) {
        e.stopPropagation();
        _filterState.text = "";
        _filterState.chartLabel = null;
        _filterState.chartViewMode = null;
        _filterState.edgeMode = "all";

        $("#bgc-filter").val("");
        $("#chart-filter-buttons .chart-filter-active").removeClass("chart-filter-active");
        $("#chart-filter-buttons .chart-filter-all").addClass("chart-filter-active");
        $("#adv-edge-toggle .adv-toggle-btn").removeClass("adv-toggle-active");
        $('#adv-edge-toggle [data-edge="all"]').addClass("adv-toggle-active");

        _filterState.similarityMin = _defaults.similarityMin;
        _filterState.similarityMax = _defaults.similarityMax;
        if (resetTypeMS) resetTypeMS();
        if (resetCatMS) resetCatMS();

        applyAllFilters();
    });

    // =========================================================================
    //  Column resize handles
    // =========================================================================
    (function () {
        var MIN_COL_WIDTH = 60;
        var $table = $("#bgc-table");
        var $ths = $table.find("thead th");

        var $resizable = $ths.filter(function () {
            return !$(this).hasClass("bgc-select-header") &&
                   !$(this).hasClass("bgc-navigate-header");
        });

        $resizable.each(function (idx) {
            var $th = $(this);
            var $handle = $('<div class="col-resize-handle"></div>');
            $th.append($handle);

            $handle.on("mousedown", function (e) {
                e.preventDefault();
                e.stopPropagation();

                var startX = e.pageX;
                var startWidth = $th.outerWidth();

                var colWidths = [];
                $resizable.each(function () { colWidths.push($(this).outerWidth()); });
                var totalBudget = 0;
                for (var k = 0; k < colWidths.length; k++) totalBudget += colWidths[k];

                var neighbourIdx = (idx < $resizable.length - 1) ? idx + 1 : idx - 1;
                var neighbourStart = colWidths[neighbourIdx];

                $(document).on("mousemove.colresize", function (ev) {
                    var delta = ev.pageX - startX;
                    var maxGrow = totalBudget - (($resizable.length - 1) * MIN_COL_WIDTH) - MIN_COL_WIDTH;
                    var newWidth = Math.max(MIN_COL_WIDTH, Math.min(startWidth + delta, startWidth + maxGrow));
                    var actualDelta = newWidth - startWidth;
                    var newNeighbour = Math.max(MIN_COL_WIDTH, neighbourStart - actualDelta);
                    actualDelta = neighbourStart - newNeighbour;
                    newWidth = startWidth + actualDelta;
                    if (newWidth < MIN_COL_WIDTH) newWidth = MIN_COL_WIDTH;

                    $th.css("width", newWidth + "px");
                    $resizable.eq(neighbourIdx).css("width", newNeighbour + "px");
                });

                $(document).on("mouseup.colresize", function () {
                    $(document).off("mousemove.colresize mouseup.colresize");
                });
            });

            $handle.on("click", function (e) {
                e.stopPropagation();
            });
        });
    })();

    // =========================================================================
    //  Custom tooltip for truncated cells
    // =========================================================================
    (function () {
        var $tooltip = null;
        var showTimer = null;

        function getTooltip() {
            if (!$tooltip) {
                $tooltip = $('<div class="dashboard-tooltip"></div>').appendTo("body").hide();
            }
            return $tooltip;
        }

        $("#bgc-table").on("mouseenter", "[data-tooltip]", function () {
            var $cell = $(this);
            var el = $cell[0];
            clearTimeout(showTimer);

            showTimer = setTimeout(function () {
                if (el.scrollWidth <= el.clientWidth) {
                    return;
                }
                var tip = getTooltip();
                tip.text($cell.attr("data-tooltip"));
                var offset = $cell.offset();
                tip.css({
                    top: offset.top + $cell.outerHeight() + 4,
                    left: offset.left
                }).show();
            }, 300);
        });

        $("#bgc-table").on("mouseleave", "[data-tooltip]", function () {
            clearTimeout(showTimer);
            if ($tooltip) {
                $tooltip.hide();
            }
        });
    })();

    // =========================================================================
    //  BGC Selection logic
    // =========================================================================
    function getSelectedAnchors() {
        var anchors = [];
        for (var key in _selectedAnchors) {
            if (_selectedAnchors[key]) {
                anchors.push(key);
            }
        }
        return anchors;
    }

    function updateSelectionUI() {
        var selected = getSelectedAnchors();
        var count = selected.length;
        var $bar = $("#bgc-selection-bar");
        var $countLabel = $("#bgc-selection-count");
        var $selectAll = $("#bgc-select-all");

        if (count > 0) {
            $bar.show();
            $countLabel.text(count + " BGC" + (count !== 1 ? "s" : "") + " selected");
        } else {
            $bar.hide();
        }

        // Manage select-all checkbox state based on current page rows
        var $visibleCheckboxes = $("#bgc-table-body .bgc-select-checkbox");
        var visibleCount = $visibleCheckboxes.length;
        var checkedVisible = $visibleCheckboxes.filter(":checked").length;

        if (visibleCount === 0 || checkedVisible === 0) {
            $selectAll.prop("checked", false);
            $selectAll.prop("indeterminate", false);
        } else if (checkedVisible === visibleCount) {
            $selectAll.prop("checked", true);
            $selectAll.prop("indeterminate", false);
        } else {
            $selectAll.prop("checked", false);
            $selectAll.prop("indeterminate", true);
        }
    }

    // Select-all checkbox: toggles all currently filtered rows (not just current page)
    $("#bgc-select-all").on("change", function () {
        var checked = $(this).prop("checked");
        for (var i = 0; i < _filteredRows.length; i++) {
            if (checked) {
                _selectedAnchors[_filteredRows[i].anchor_id] = true;
            } else {
                delete _selectedAnchors[_filteredRows[i].anchor_id];
            }
        }
        renderPage();
    });

    // Individual checkbox change (delegated)
    $("#bgc-table-body").on("change", ".bgc-select-checkbox", function () {
        var anchor = $(this).data("anchor");
        if ($(this).prop("checked")) {
            _selectedAnchors[anchor] = true;
        } else {
            delete _selectedAnchors[anchor];
        }
        updateSelectionUI();
    });

    // "View selected" button
    $("#bgc-view-selected").on("click", function () {
        var anchors = [];
        for (var i = 0; i < _filteredRows.length; i++) {
            if (_selectedAnchors[_filteredRows[i].anchor_id]) {
                anchors.push(_filteredRows[i].anchor_id);
            }
        }
        if (anchors.length === 0) {
            return;
        }
        sessionStorage.setItem("antismash_bgc_selection", JSON.stringify(anchors));
        window.location.href = "regions.html";
    });

    // "Clear selection" button
    $("#bgc-clear-selection").on("click", function () {
        _selectedAnchors = {};
        $("#bgc-select-all").prop("checked", false).prop("indeterminate", false);
        renderPage();
    });

    // =========================================================================
    //  Initial render
    // =========================================================================
    renderPage();
});
