/* Lazy-loading module for per-region data and HTML fragments.
 *
 * Region data is stored in compressed bundles (regions_data/0001.js, ...) of a
 * few MB each, which add anchor -> base64-encoded gzipped JSON entries to
 * _regionsData. regions_data/index.js maps each anchor to its bundle. Only the
 * bundle of a region being shown is loaded (via <script>, so file:// works) and
 * only that region is decompressed, using DecompressionStream (with pako fallback).
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

    // Anchors whose currently injected HTML has had viewer handlers bound
    var _bound = {};

    // Bundle loading state: number -> "loading" | "loaded", plus waiting callbacks
    var _bundles = {};
    var _bundleCallbacks = {};

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
        // fresh DOM, so any previously bound handlers are gone
        _bound[anchor] = false;
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

    /** Create the hidden container a region's HTML is injected into, unless it exists.
     *  Returns whether it was created. */
    function _ensureContainer(anchor) {
        if (document.getElementById(anchor)) {
            return false;
        }
        var div = document.createElement("div");
        div.className = "page";
        div.id = anchor;
        div.setAttribute("data-lazy", "true");
        div.style.display = "none";
        div.innerHTML = '<div class="lazy-placeholder"><div class="loading-spinner"></div>'
            + '<p>Loading region ' + anchor.replace(/^r(\d+)c(\d+)$/, "$1.$2") + '...</p></div>';
        var marker = document.getElementById("lazy-placeholder-anchor");
        if (marker) {
            marker.parentNode.insertBefore(div, marker.nextSibling);
        } else {
            document.body.appendChild(div);
        }
        return true;
    }

    /** Load a regions_data bundle once, then call the callback (ok: true/false) */
    function _loadBundle(bundle, callback) {
        if (_bundles[bundle] === "loaded") {
            callback(true);
            return;
        }
        _bundleCallbacks[bundle] = _bundleCallbacks[bundle] || [];
        _bundleCallbacks[bundle].push(callback);
        if (_bundles[bundle] === "loading") {
            return;
        }
        _bundles[bundle] = "loading";
        var finish = function (ok) {
            _bundles[bundle] = ok ? "loaded" : undefined;
            var cbs = _bundleCallbacks[bundle];
            _bundleCallbacks[bundle] = [];
            for (var i = 0; i < cbs.length; i++) {
                cbs[i](ok);
            }
        };
        var script = document.createElement("script");
        script.src = "regions_data/" + ("000" + bundle).slice(-Math.max(4, String(bundle).length)) + ".js";
        script.onload = function () { finish(true); };
        script.onerror = function () { finish(false); };
        document.head.appendChild(script);
    }

    /** Core loading routine for a single anchor */
    function _load(anchor, callback) {
        if (_ensureContainer(anchor) && _state[anchor] === "loaded") {
            // its container had been removed (detailed pagination): inject it again
            _injectHTML(anchor);
        }
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

        if (typeof _regionsData === "undefined" || !_regionsData[anchor]) {
            var bundle = typeof _regionsDataIndex !== "undefined" ? _regionsDataIndex[anchor] : undefined;
            if (bundle === undefined) {
                console.warn("No data found for region: " + anchor);
                _state[anchor] = "pending";
                return;
            }
            _loadBundle(bundle, function (ok) {
                if (!ok || typeof _regionsData === "undefined" || !_regionsData[anchor]) {
                    console.warn("Failed to load data for region: " + anchor);
                    _state[anchor] = "pending";
                    return;
                }
                _decode(anchor);
            });
            return;
        }
        _decode(anchor);
    }

    /** Decompress a loaded region's data and inject it into the page */
    function _decode(anchor) {
        var bytes = _b64ToBytes(_regionsData[anchor]);

        _gunzip(bytes).then(function (decompressed) {
            var jsonStr = new TextDecoder().decode(decompressed);
            var data = JSON.parse(jsonStr);

            // Populate the global data structures; regions without module
            // results still need an entry, which the viewer expects
            _regions[anchor] = data.region;
            _results[anchor] = data.results || {};
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
    }

    /** Bind the viewer's handlers for one region's freshly injected HTML.
     */
    function _bindRegion(anchor) {
        var container = document.getElementById(anchor);
        var real$ = window.$;
        var scoped$ = function (selector, context) {
            if (typeof selector === "string" && context === undefined && selector.charAt(0) !== "<") {
                // search only the region's subtree: a page-wide query is slow on large
                // pages. Native querySelectorAll (unlike jQuery's find) also matches
                // selectors like "#<anchor> * .x" that start at the container itself.
                try {
                    var nodes = Array.prototype.slice.call(container.querySelectorAll(selector));
                    if (container.matches(selector)) {
                        nodes.unshift(container);
                    }
                    return real$(nodes);
                } catch (e) {
                    // not plain CSS (jQuery-only syntax): filter a page-wide query
                    return real$(selector).filter(function () {
                        return this === container || container.contains(this);
                    });
                }
            }
            return real$(selector, context);
        };
        real$.extend(scoped$, real$);
        var regionOnly = Object.create(_regions);
        regionOnly.order = [anchor];
        window.$ = scoped$;
        try {
            viewer["start"](regionOnly, _results, []);
        } finally {
            window.$ = real$;
            // restore the full region set (next/prev navigation uses its order)
            viewer["setData"](_regions, _results);
        }
        _bound[anchor] = true;
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
        if (typeof viewer === "undefined" || !viewer["start"]) {
            return;
        }
        if (!_bound[anchor] && document.getElementById(anchor)) {
            _bindRegion(anchor);
        } else {
            viewer["switchToRegion"]();
        }
    }

    /** Step to the previous/next region in _regions.order, as the viewer's arrows do:
     *  from the overview to the first/last region, and past either end to the overview */
    function _step(direction) {
        var order = _regions.order;
        var current = window.location.hash.replace(/^#/, "");
        var target;
        if (!current || current === "overview") {
            target = direction > 0 ? order[0] : order[order.length - 1];
        } else {
            target = order[order.indexOf(current) + direction];
        }
        if (target === undefined) {
            // the lazy loader ignores an empty hash, so show the overview directly
            window.location.hash = "";
            viewer["switchToRegion"]();
            return;
        }
        window.location.hash = target;
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

        // Prev/next arrows
        $("#prev-region, #next-region").on("click", function (event) {
            event.stopImmediatePropagation();
            _step(this.id === "next-region" ? 1 : -1);
        });
        window.addEventListener("keyup", function (event) {
            var key = event.keyCode;
            var direction = (key === 37 || key === 87) ? -1 : ((key === 39 || key === 69) ? 1 : 0);
            if (!direction) return;
            event.stopPropagation();
            _step(direction);
        }, true);

        // Hook into overview table row clicks
        $(document).on("click", ".linked-row", function () {
            var dataAnchor = $(this).data("anchor") || "";
            var a = dataAnchor.replace(/^#/, "");
            if (a && _state[a] !== undefined) {
                loadAndShow(a);
            }
        });

        // Bind the page-level viewer handlers. No records are passed, as the
        // overview minimaps are drawn for the visible overview page only.
        if (typeof viewer !== "undefined" && viewer["start"]) {
            viewer["start"](_regions, _results, []);
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
