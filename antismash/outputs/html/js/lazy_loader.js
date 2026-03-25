/* Lazy-loading module for per-region data and HTML fragments.
 *
 * Region data is stored in a single compressed file (regions_data.js) as a
 * JS object mapping anchor -> base64-encoded gzipped JSON. Each region is
 * decompressed on demand using DecompressionStream (with pako fallback).
 *
 * This replaces the previous approach of one <script> per region, reducing
 * thousands of tiny files to a single data file while keeping file:// support.
 */
var lazyRegionLoader = (function () {
    "use strict";

    var _regions = null;   // reference to global all_regions
    var _records = null;   // reference to global recordData
    var _results = null;   // reference to global resultsData

    // Track loading state per anchor: "pending" | "loading" | "loaded"
    var _state = {};

    // Callbacks waiting for a specific anchor to finish loading
    var _callbacks = {};

    // O(1) anchor -> {recordIndex, regionIndex} map, built once in init()
    var _anchorMap = {};

    // Whether the browser supports DecompressionStream
    var _hasDecompressionStream = typeof DecompressionStream !== "undefined";

    /** Decode a base64 string to a Uint8Array */
    function _b64ToBytes(b64) {
        var binary = atob(b64);
        var bytes = new Uint8Array(binary.length);
        for (var i = 0; i < binary.length; i++) {
            bytes[i] = binary.charCodeAt(i);
        }
        return bytes;
    }

    /** Decompress gzipped bytes using DecompressionStream API */
    function _decompressStream(bytes) {
        var ds = new DecompressionStream("gzip");
        var writer = ds.writable.getWriter();
        writer.write(bytes);
        writer.close();
        var reader = ds.readable.getReader();
        var chunks = [];
        function pump() {
            return reader.read().then(function (result) {
                if (result.done) {
                    var totalLen = 0;
                    for (var i = 0; i < chunks.length; i++) totalLen += chunks[i].length;
                    var combined = new Uint8Array(totalLen);
                    var offset = 0;
                    for (var j = 0; j < chunks.length; j++) {
                        combined.set(chunks[j], offset);
                        offset += chunks[j].length;
                    }
                    return combined;
                }
                chunks.push(result.value);
                return pump();
            });
        }
        return pump();
    }

    /** Decompress gzipped bytes using pako fallback */
    function _decompressPako(bytes) {
        return Promise.resolve(pako.inflate(bytes));
    }

    /** Decompress gzipped bytes, choosing the best available method */
    function _gunzip(bytes) {
        if (_hasDecompressionStream) {
            return _decompressStream(bytes);
        }
        if (typeof pako !== "undefined") {
            return _decompressPako(bytes);
        }
        return Promise.reject(new Error("No decompression method available"));
    }

    /** Build the anchor lookup map from recordData */
    function _buildAnchorMap() {
        _anchorMap = {};
        for (var i = 0; i < _records.length; i++) {
            for (var j = 0; j < _records[i].regions.length; j++) {
                _anchorMap[_records[i].regions[j].anchor] = {recordIndex: i, regionIndex: j};
            }
        }
    }

    /** Find the record entry in recordData that contains a region with the given anchor */
    function _findRecordForAnchor(anchor) {
        return _anchorMap[anchor] || null;
    }

    /** Patch the lightweight region in recordData with the full data after load */
    function _patchRecordData(anchor) {
        var info = _findRecordForAnchor(anchor);
        if (info && _regions[anchor]) {
            _records[info.recordIndex].regions[info.regionIndex] = _regions[anchor];
        }
    }

    /** Inject the HTML content from regionHTML into the placeholder div */
    function _injectHTML(anchor) {
        var container = document.getElementById(anchor);
        if (!container) {
            return;
        }
        if (typeof regionHTML !== "undefined" && regionHTML[anchor]) {
            container.innerHTML = regionHTML[anchor];
        } else {
            // Minimal fallback — just provide an SVG container
            container.innerHTML =
                '<div class="region-grid"><div class="content">' +
                '<div class="description-container">' +
                '<div class="heading">Region ' + anchor + '</div>' +
                '<div class="region-svg-container"><div id="' + anchor + '-svg"></div></div>' +
                '</div></div>' +
                '<div class="focus-panel"><div class="heading"><span>Gene details</span></div>' +
                '<div class="focus-panel-content focus-panel-content-' + anchor + '">' +
                '<div style="text-align:center;margin-top:30%;">Select a gene to view details</div>' +
                '</div></div></div>';
        }
    }

    /** Fire all waiting callbacks for an anchor */
    function _fireCallbacks(anchor) {
        var cbs = _callbacks[anchor] || [];
        _callbacks[anchor] = [];
        for (var i = 0; i < cbs.length; i++) {
            cbs[i]();
        }
    }

    /** Core loading routine for a single anchor */
    function _load(anchor, callback) {
        if (_state[anchor] === "loaded") {
            if (callback) { callback(); }
            return;
        }
        if (_state[anchor] === "loading") {
            if (callback) {
                _callbacks[anchor] = _callbacks[anchor] || [];
                _callbacks[anchor].push(callback);
            }
            return;
        }

        _state[anchor] = "loading";
        if (callback) {
            _callbacks[anchor] = [callback];
        }

        // Check if compressed data is available
        if (typeof _regionsData !== "undefined" && _regionsData[anchor]) {
            var b64 = _regionsData[anchor];
            var bytes = _b64ToBytes(b64);

            _gunzip(bytes).then(function (decompressed) {
                var jsonStr = new TextDecoder().decode(decompressed);
                var data = JSON.parse(jsonStr);

                // Populate the global data structures
                _regions[anchor] = data.region;
                if (data.results) {
                    _results[anchor] = data.results;
                }
                if (data.html) {
                    regionHTML[anchor] = data.html;
                }

                // Free the base64 string to reduce memory
                delete _regionsData[anchor];

                _patchRecordData(anchor);
                _injectHTML(anchor);
                _state[anchor] = "loaded";
                _fireCallbacks(anchor);
            }).catch(function (err) {
                console.warn("Failed to decompress data for region: " + anchor, err);
                _state[anchor] = "pending";
            });
        } else {
            // Fallback: try loading individual file via script injection
            var script = document.createElement("script");
            script.src = "regions/" + anchor + ".js";
            script.onload = function () {
                _patchRecordData(anchor);
                _injectHTML(anchor);
                _state[anchor] = "loaded";
                _fireCallbacks(anchor);
            };
            script.onerror = function () {
                console.warn("Failed to load data for region: " + anchor);
                _state[anchor] = "pending";
            };
            document.head.appendChild(script);
        }
    }

    /** Draw the SVG for a single region using the viewer (debounced) */
    var _lastDrawnAnchor = null;
    var _lastDrawnTime = 0;

    function _drawRegion(anchor) {
        var now = Date.now();
        if (anchor === _lastDrawnAnchor && now - _lastDrawnTime < 200) {
            return;
        }
        _lastDrawnAnchor = anchor;
        _lastDrawnTime = now;
        if (typeof viewer !== "undefined" && viewer["start"]) {
            viewer["start"](_regions, _results, _records);
        }
    }

    /** Navigate to (show) a specific region, loading it first if needed */
    function loadAndShow(anchor) {
        _load(anchor, function () {
            _drawRegion(anchor);
        });
    }

    /** Initialise the lazy loader, hook into hash navigation */
    function init(regions, records, results) {
        _regions = regions;
        _records = records;
        _results = results;

        // Build O(1) anchor lookup map
        _buildAnchorMap();

        // Mark all known anchors as pending
        for (var i = 0; i < _regions.order.length; i++) {
            _state[_regions.order[i]] = "pending";
        }

        // If the URL already has a hash pointing to a region, load it
        var hash = window.location.hash.replace(/^#/, "");
        if (hash && _state[hash] !== undefined) {
            loadAndShow(hash);
        }

        // Listen for hash changes to trigger lazy loading
        $(window).on("hashchange", function () {
            var h = window.location.hash.replace(/^#/, "");
            if (h && _state[h] !== undefined) {
                loadAndShow(h);
            }
        });

        // Hook into region button clicks
        $("#regionbuttons").on("click", "a", function () {
            var href = $(this).attr("href") || "";
            var a = href.replace(/^#/, "");
            if (a && _state[a] !== undefined) {
                loadAndShow(a);
            }
        });

        // Hook into prev/next navigation arrows
        $("#prev-region, #next-region").on("click", function () {
            // Let antismash.js update the hash first, then check
            setTimeout(function () {
                var h = window.location.hash.replace(/^#/, "");
                if (h && _state[h] !== undefined) {
                    loadAndShow(h);
                }
            }, 50);
        });

        // Hook into overview table row clicks
        $(document).on("click", ".linked-row", function () {
            var dataAnchor = $(this).data("anchor") || "";
            var a = dataAnchor.replace(/^#/, "");
            if (a && _state[a] !== undefined) {
                loadAndShow(a);
            }
        });

        // Call viewer start with the lightweight data for the overview
        if (typeof viewer !== "undefined" && viewer["start"]) {
            viewer["start"](_regions, _results, _records);
        }
    }

    // Preload queue state
    var _preloadCancelled = false;

    /** Cancel any in-flight preload queue (safe to call even if none running) */
    function cancelPreload() {
        _preloadCancelled = true;
    }

    /**
     * Pre-load a batch of anchors with a concurrency limit.
     * @param {string[]} anchors - anchors to preload
     * @param {number} [concurrency=6] - max parallel loads
     */
    function preloadAnchors(anchors, concurrency) {
        if (!anchors || anchors.length === 0) return;
        concurrency = concurrency || 6;
        _preloadCancelled = false;

        var queue = [];
        for (var i = 0; i < anchors.length; i++) {
            if (_state[anchors[i]] === "pending") {
                queue.push(anchors[i]);
            }
        }

        var idx = 0;
        var active = 0;

        function next() {
            while (active < concurrency && idx < queue.length) {
                if (_preloadCancelled) return;
                (function (anchor) {
                    active++;
                    _load(anchor, function () {
                        active--;
                        if (!_preloadCancelled) next();
                    });
                })(queue[idx++]);
            }
        }
        next();
    }

    /** Check if a region's data has been loaded */
    function isLoaded(anchor) {
        return _state[anchor] === "loaded";
    }

    /** Re-inject cached HTML for already-loaded anchors whose DOM was rebuilt */
    function reinjectLoaded(anchors) {
        for (var i = 0; i < anchors.length; i++) {
            var anchor = anchors[i];
            if (_state[anchor] === "loaded") {
                _injectHTML(anchor);
            }
        }
    }

    /** Re-scan _regions.order and register any new anchors not yet tracked */
    function _reinitState() {
        if (!_regions) return;
        for (var i = 0; i < _regions.order.length; i++) {
            var anchor = _regions.order[i];
            if (!_state[anchor]) {
                _state[anchor] = "pending";
            }
        }
    }

    return {
        init: init,
        loadAndShow: loadAndShow,
        preloadAnchors: preloadAnchors,
        cancelPreload: cancelPreload,
        isLoaded: isLoaded,
        reinjectLoaded: reinjectLoaded,
        _reinitState: _reinitState
    };
})();
