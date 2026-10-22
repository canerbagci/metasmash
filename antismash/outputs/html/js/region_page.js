/* Standalone single-region page (region.html#<anchor>), opened from the dashboard.
 *
 * Loads only the viewer and the regions_data bundle holding that region (via the
 * lazy loader), so opening a BGC doesn't pay for the full regions.html overview.
 * Prev/next step through the dashboard's current filtered + sorted order, which
 * the dashboard stores in sessionStorage when a region link is clicked.
 */
(function () {
    "use strict";

    var NAV_KEY = "metasmash_region_nav";
    var ANCHOR_PATTERN = /^r\d+c\d+$/;

    // each region gets a fresh page, so the viewer is only ever started once
    window.addEventListener("hashchange", function () {
        window.location.reload();
    });

    var anchor = window.location.hash.replace(/^#/, "");
    // the viewer navigates past the first/last region to "#overview"
    if (anchor === "" || anchor === "overview") {
        window.location.replace("index.html");
        return;
    }

    function showMessage(text) {
        $("#region-page-message").text(text).show();
    }

    function label(anchorId) {
        var match = anchorId.match(/^r(\d+)c(\d+)$/);
        return match ? match[1] + "." + match[2] : anchorId;
    }

    if (!ANCHOR_PATTERN.test(anchor)) {
        $("#region-page-nav").hide();
        showMessage("Unknown region: " + anchor);
        return;
    }

    $("#region-page-detailed-link").attr("href", "regions.html#" + anchor);

    var nav = null;
    try {
        nav = JSON.parse(sessionStorage.getItem(NAV_KEY) || "null");
    } catch (e) {
        nav = null;
    }
    var position = Array.isArray(nav) ? nav.indexOf(anchor) : -1;
    if (position === -1) {
        nav = null;
        $("#region-page-info").text("Region " + label(anchor));
    } else {
        $("#region-page-info").text("Region " + label(anchor) + " (" + (position + 1)
                                    + " of " + nav.length + ")");
        var linkTo = function (selector, target) {
            if (!target) return;
            $(selector).prop("disabled", false).on("click", function () {
                window.location.hash = target;
            });
        };
        linkTo("#region-page-prev", nav[position - 1]);
        linkTo("#region-page-next", nav[position + 1]);
    }

    // The viewer's arrow-key navigation (left/W, right/E) would try to draw the next
    // region before it is loaded, so handle those keys here, before they reach it.
    window.addEventListener("keyup", function (event) {
        var key = event.keyCode;
        var step = (key === 37 || key === 87) ? -1 : ((key === 39 || key === 69) ? 1 : 0);
        if (!step) return;
        // also keeps the lazy loader's own arrow-key navigation from running
        event.stopImmediatePropagation();
        if (nav && nav[position + step]) {
            window.location.hash = nav[position + step];
        }
    }, true);

    if (typeof _regionsDataIndex === "undefined" || _regionsDataIndex[anchor] === undefined) {
        showMessage("No data found for region " + label(anchor) + ".");
        return;
    }

    // placeholder the lazy loader injects the region's HTML into
    var page = document.createElement("div");
    page.className = "page";
    page.id = anchor;
    document.getElementById("region-page-container").appendChild(page);

    // the viewer expects a results entry for the region shown, even before it is loaded
    resultsData[anchor] = {};
    all_regions.order = [anchor];
    // loads the region named in the URL hash, then binds and draws it
    lazyRegionLoader.init(all_regions, [], resultsData);
})();
