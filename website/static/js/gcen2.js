/* GCEN2 website: colour mode switch, copy buttons, gallery filter and dialog.
   The colour mode is first set by the inline script in <head> (no flash). */
(function () {
  "use strict";

  var KEY = "gcen2-theme";

  function stored() {
    try { return localStorage.getItem(KEY); } catch (e) { return null; }
  }
  function store(v) {
    try { localStorage.setItem(KEY, v); } catch (e) { /* storage blocked: keep for this page */ }
  }
  function preferred() {
    var s = stored();
    if (s === "light" || s === "dark") return s;
    return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }
  function apply(mode) {
    document.documentElement.setAttribute("data-bs-theme", mode);
    // the introduction video in the colours of the mode, unless it is already playing
    document.querySelectorAll("video[data-gc-video-light]").forEach(function (v) {
      if (!v.paused || v.currentTime > 0) return;
      var which = mode === "dark" ? "dark" : "light";
      var src = v.getAttribute("data-gc-video-" + which);
      if (v.getAttribute("src") !== src) {
        v.setAttribute("src", src);
        v.setAttribute("poster", v.getAttribute("data-gc-poster-" + which));
      }
    });
    document.querySelectorAll("[data-gc-theme-toggle]").forEach(function (b) {
      var dark = mode === "dark";
      b.setAttribute("aria-pressed", dark ? "true" : "false");
      b.setAttribute("aria-label", dark ? "Switch to light mode" : "Switch to dark mode");
      b.setAttribute("title", dark ? "Light mode" : "Dark mode");
      b.querySelectorAll("[data-gc-icon]").forEach(function (i) {
        i.classList.toggle("d-none", i.getAttribute("data-gc-icon") !== (dark ? "sun" : "moon"));
      });
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    apply(preferred());
    document.querySelectorAll("[data-gc-theme-toggle]").forEach(function (b) {
      b.addEventListener("click", function () {
        var next = document.documentElement.getAttribute("data-bs-theme") === "dark" ? "light" : "dark";
        store(next);
        apply(next);
      });
    });
    window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", function () {
      if (!stored()) apply(preferred());
    });

    // The play button over the video poster: until the video first plays, it replaces
    // the browser's controls (without JavaScript the controls stay, and the button hidden).
    document.querySelectorAll("[data-gc-video-play]").forEach(function (btn) {
      var v = btn.parentElement.querySelector("video");
      v.controls = false;
      btn.hidden = false;
      btn.addEventListener("click", function () { v.controls = true; v.play(); });
      v.addEventListener("play", function () { v.controls = true; btn.hidden = true; });
    });

    // Copy buttons: copy the text of the <pre> next to the button, without prompts.
    document.querySelectorAll(".gc-copy").forEach(function (btn) {
      btn.addEventListener("click", function () {
        var pre = btn.parentElement.querySelector("pre");
        var clone = pre.cloneNode(true);
        clone.querySelectorAll(".gc-prompt").forEach(function (p) { p.remove(); });
        var text = clone.textContent.replace(/\s+$/, "") + "\n";
        var done = function () {
          var label = btn.querySelector(".gc-copy-label");
          btn.classList.add("copied");
          if (label) label.textContent = "Copied";
          setTimeout(function () {
            btn.classList.remove("copied");
            if (label) label.textContent = "Copy";
          }, 1500);
        };
        var fallback = function () {
          var r = document.createRange();
          r.selectNodeContents(pre);
          var s = window.getSelection();
          s.removeAllRanges();
          s.addRange(r);
          try { document.execCommand("copy"); done(); } catch (e) { /* text stays selected */ }
        };
        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(text).then(done, fallback);
        } else {
          fallback();
        }
      });
    });

    // Gallery filter
    var filter = document.querySelector("[data-gc-filter]");
    if (filter) {
      var items = document.querySelectorAll("[data-gc-kind]");
      var status = document.getElementById("galleryStatus");
      filter.querySelectorAll("button").forEach(function (b) {
        b.addEventListener("click", function () {
          var kind = b.getAttribute("data-gc-show");
          filter.querySelectorAll("button").forEach(function (o) {
            var on = o === b;
            o.classList.toggle("active", on);
            o.setAttribute("aria-pressed", on ? "true" : "false");
          });
          var n = 0;
          items.forEach(function (it) {
            var show = kind === "all" || it.getAttribute("data-gc-kind") === kind;
            it.classList.toggle("d-none", !show);
            if (show) n += 1;
          });
          if (status) status.textContent = n + (n === 1 ? " figure" : " figures");
        });
      });
    }

    // Gallery dialog
    var modal = document.getElementById("figureModal");
    if (modal) {
      modal.addEventListener("show.bs.modal", function (ev) {
        var t = ev.relatedTarget;
        modal.querySelector(".modal-title").textContent = t.getAttribute("data-gc-title");
        var img = modal.querySelector(".modal-body img");
        img.src = t.getAttribute("data-gc-src");
        img.alt = t.getAttribute("data-gc-alt");
        img.width = +t.getAttribute("data-gc-w");
        img.height = +t.getAttribute("data-gc-h");
        modal.querySelector("[data-gc-modal-cmd]").textContent = t.getAttribute("data-gc-cmd");
      });
    }
  });
})();
