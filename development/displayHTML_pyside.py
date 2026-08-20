"""
displayHTML_pyside.py
──────────────────────
Loads a simpleNMR D3 HTML file in a QWebEngineView and intercepts the
"Export" button so the JSON data is sent back to Python via QWebChannel
instead of being saved to disk.

No modifications to the HTML file are required.

JS console.log/warn/error calls from the page (including the injected
script below) are forwarded to Python's stdout via a QWebEnginePage
subclass that overrides javaScriptConsoleMessage() — without this
override, JS console output is invisible from PySide6/Qt.

Usage:
    python displayHTML_pyside.py path/to/Cytochalasin-B-kate_d3.html

Requirements:
    pip install PySide6
    (PySide6 bundles QtWebEngine — no separate WebEngine package needed)
"""

import sys
import os
import json
from datetime import datetime

# ── QWebEngineWidgets MUST be imported before QApplication ──────────────────
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebEngineCore import QWebEngineSettings, QWebEnginePage
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWidgets import QApplication, QMainWindow, QStatusBar, QMessageBox
from PySide6.QtCore import QObject, QUrl, Signal, Slot
from PySide6.QtGui import QDesktopServices

# The clipboard-permission API changed between PySide6 releases (older:
# QWebEnginePage.featurePermissionRequested / Feature enum; newer:
# QWebEnginePage.permissionRequested / QWebEnginePermission object). We
# probe for whichever is present at import time rather than assuming one,
# so this file works across the version you had before and PySide6 on
# Python 3.13 / arm64.
try:
    from PySide6.QtWebEngineCore import QWebEnginePermission  # newer API
    HAS_NEW_PERMISSION_API = True
except ImportError:
    HAS_NEW_PERMISSION_API = False

# Helps prevent sandbox errors on macOS.
#
# --single-process: added after diagnostics showed the crash only happens
# when this script is launched by JASON, not from Terminal. The key
# difference is XPC_SERVICE_NAME - JASON-launched processes run inside
# JASON's own app-scoped Mach/XPC bootstrap namespace rather than the
# normal session-wide one a Terminal shell gets. Chromium's multi-process
# architecture (browser/GPU/renderer processes) hands off Mach ports
# between those processes - exactly what "mach_msg receive: (ipc/rcv) msg
# too large" in mojo/core/channel_mac.cc is failing to do. Running
# QtWebEngine single-process removes that cross-process handoff entirely,
# so it shouldn't matter which Mach bootstrap namespace we're in.
#
# Trade-off: single-process mode gives up Chromium's process isolation -
# a renderer crash now takes the whole app down instead of being caught by
# on_render_process_terminated() below. For this single-purpose viewer
# window, that trade-off is acceptable.
#
# (The earlier --disable-gpu-compositing flag has been removed: the
# diagnostics comparison showed identical code/content crashing under
# JASON but working from Terminal, which points at the process-launch
# environment rather than GPU/content, so that flag wasn't addressing the
# actual cause.)
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--no-sandbox --single-process"

# ── JavaScript injected after the page has finished loading ─────────────────
# Strategy:
#   1. Dynamically load qwebchannel.js from Qt's built-in resource path.
#   2. Once loaded, open the channel and store window.pyBridge.
#   3. Patch HTMLAnchorElement.prototype.click so that any programmatic
#      anchor.click() on a data:text/json URL is intercepted and the JSON
#      is sent to Python instead of triggering a browser download.
#      If pyBridge is not available (standalone browser use), the original
#      click behaviour is preserved as a fallback.
INJECT_JS = r"""
(function () {
    try {
        // ── 1. Load qwebchannel.js from Qt's built-in resource ──────────
        var script = document.createElement('script');
        script.src = 'qrc:///qtwebchannel/qwebchannel.js';

        script.onload = function () {
            // ── 2. Open the channel Qt registered as "pyBridge" ─────────
            try {
                new QWebChannel(qt.webChannelTransport, function (channel) {
                    window.pyBridge = channel.objects.pyBridge;
                    if (!window.pyBridge) {
                        console.error('[nmr_viewer] QWebChannel opened but pyBridge is undefined!');
                    }
                });
            } catch (e) {
                console.error('[nmr_viewer] Failed to open QWebChannel: ' + e);
            }
        };

        script.onerror = function () {
            console.error('[nmr_viewer] Failed to load qwebchannel.js from qrc:///qtwebchannel/qwebchannel.js');
        };

        document.head.appendChild(script);

        // ── 3. Patch anchor clicks to intercept JSON data-URI downloads ──
        //    We patch the prototype once, before any user interaction occurs.
        //    exportToMnova() calls anchor.click() programmatically, so this
        //    fires reliably IF the export mechanism uses <a>.click().
        var _origClick = HTMLAnchorElement.prototype.click;

        HTMLAnchorElement.prototype.click = function () {
            var href = this.href || '';

            if (href.startsWith('data:text/json')) {
                if (window.pyBridge) {
                    // Decode the JSON from the data URI and send to Python
                    // — this calls DataBridge.receiveExportData().
                    var encoded = href.replace(/^data:text\/json;charset=utf-8,/, '');
                    var jsonStr = decodeURIComponent(encoded);
                    window.pyBridge.receiveExportData(jsonStr);

                    // Do NOT call _origClick — skip the file download
                    return;
                } else {
                    console.error('[nmr_viewer] href matched data:text/json but pyBridge is NOT available yet. ' +
                                   'Falling back to normal download. (Was the page loaded before the channel opened?)');
                }
            }

            // Anything else (non-JSON links, pyBridge not ready, etc.) behaves normally
            _origClick.call(this);
        };

        // ── 3b. Guard against navigator.clipboard.writeText() hanging ────
        //    exportToMnova() does:
        //        await navigator.clipboard.writeText(workingDirectory);
        //    before it ever builds the export <a> and calls .click().
        //    QtWebEngine frequently never resolves OR rejects this promise
        //    (no permission prompt implementation), which means the await
        //    blocks forever and the anchor/.click() below it never runs —
        //    this looks exactly like "the button click isn't intercepted"
        //    from the outside, when really the export code never got that
        //    far. We patch writeText so it always settles within 800ms,
        //    trying the real clipboard write in the background but never
        //    letting it block the caller.
        if (navigator.clipboard && navigator.clipboard.writeText) {
            var _origWriteText = navigator.clipboard.writeText.bind(navigator.clipboard);
            navigator.clipboard.writeText = function (text) {
                return Promise.race([
                    _origWriteText(text)
                        .catch(function (e) {
                            console.warn('[nmr_viewer] clipboard.writeText rejected (non-fatal): ' + e);
                        }),
                    new Promise(function (resolve) {
                        setTimeout(function () {
                            console.warn('[nmr_viewer] clipboard.writeText did not settle within 800ms — continuing anyway so export is not blocked.');
                            resolve();
                        }, 800);
                    }),
                ]);
            };
        } else {
            console.warn('[nmr_viewer] navigator.clipboard.writeText not available in this environment.');
        }
    } catch (e) {
        console.error('[nmr_viewer] Injection failed with exception: ' + e);
    }
})();
"""


# ── QWebEnginePage subclass that forwards JS console output to Python ───────
class LoggingWebEnginePage(QWebEnginePage):
    """
    By default, console.log/warn/error calls made by JavaScript running in
    a QWebEngineView are NOT printed anywhere Python can see. Overriding
    javaScriptConsoleMessage() forwards them to stdout, which is essential
    for actually seeing the [nmr_viewer] trace messages from INJECT_JS.
    """

    LEVEL_NAMES = {
        QWebEnginePage.JavaScriptConsoleMessageLevel.InfoMessageLevel: "INFO",
        QWebEnginePage.JavaScriptConsoleMessageLevel.WarningMessageLevel: "WARN",
        QWebEnginePage.JavaScriptConsoleMessageLevel.ErrorMessageLevel: "ERROR",
    }

    def javaScriptConsoleMessage(self, level, message, lineNumber, sourceID):
        level_name = self.LEVEL_NAMES.get(level, str(level))
        short_source = os.path.basename(sourceID) if sourceID else "?"
        print(f"[JS:{level_name}] {message}  ({short_source}:{lineNumber})")

    def createWindow(self, _window_type):
        """Handle JS window.open() calls (e.g. the HTML's Help button).

        QWebEnginePage does nothing with window.open() unless this is
        overridden - by default the call just silently no-ops. The HTML
        template's Help button does:
            window.open('http://simplenmr.pythonanywhere.com/documentation/...', '_blank')
        which is meant to open the docs in a real browser, not a second
        embedded Qt window, so we hand the URL off to the system's default
        browser via QDesktopServices and return a disposable page for
        Chromium's window.open() machinery to write its navigation into
        (it needs *some* QWebEnginePage back, or the JS call errors).
        """
        proxy_page = QWebEnginePage(self)
        proxy_page.urlChanged.connect(
            lambda url, page=proxy_page: (QDesktopServices.openUrl(url), page.deleteLater())
        )
        return proxy_page


# ── Renderer-crash handling ──────────────────────────────────────────────────
# QtWebEngine runs each page in a separate renderer process; if that process
# crashes (e.g. the macOS "mach_msg ... msg too large" IPC failure), Qt does
# NOT raise a Python exception - the page just goes blank. renderProcessTerminated
# is the only signal that surfaces this, so without connecting it, a crash is
# silently invisible from the Python side and only visible as raw stderr
# (which is how JASON ends up showing the mojo/channel_mac.cc error text).
TERMINATION_STATUS_NAMES = {
    QWebEnginePage.RenderProcessTerminationStatus.NormalTerminationStatus: "Normal",
    QWebEnginePage.RenderProcessTerminationStatus.AbnormalTerminationStatus: "Abnormal",
    QWebEnginePage.RenderProcessTerminationStatus.CrashedTerminationStatus: "Crashed",
    QWebEnginePage.RenderProcessTerminationStatus.KilledTerminationStatus: "Killed",
}


def _log_viewer_crash(html_path: str, status_name: str, exit_code: int):
    """Append renderer-crash details to a log file next to the HTML output.

    Mirrors _log_submission_error()'s pattern in simpleNMRjeolTools_v6.py,
    but for viewer-side (renderer process) failures rather than server
    submission failures - these are different failure points and were
    previously not logged anywhere. Never raises.
    """
    try:
        log_dir = os.path.dirname(html_path)
        log_path = os.path.join(log_dir, "simpleNMR_viewer_crash.log")
        timestamp = datetime.now().isoformat(timespec="seconds")
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(f"\n{'=' * 70}\n")
            f.write(f"Timestamp: {timestamp}\n")
            f.write(f"HTML file: {html_path}\n")
            f.write(f"Termination status: {status_name}\n")
            f.write(f"Exit code: {exit_code}\n")
        return log_path
    except Exception as e:
        print(f"Failed to write viewer crash log: {e}")
        return None


# ── Python-side bridge object exposed to JavaScript ─────────────────────────
class DataBridge(QObject):
    """
    Registered with QWebChannel under the name 'pyBridge'.
    JavaScript calls  window.pyBridge.receiveExportData(jsonString)
    which triggers the receiveExportData slot below.
    """

    # Emitted with the parsed dict whenever export data arrives
    exportDataReceived = Signal(dict)

    @Slot(str)
    def receiveExportData(self, json_str: str):
        """Called from JavaScript when the Export button is pressed."""
        try:
            data = json.loads(json_str)
        except json.JSONDecodeError as e:
            print(f"[Python][DataBridge] JSON decode error: {e}")
            return

        working_fn = data.get("workingFilename", "unknown")
        n_nodes = len(data.get("nodes_now", []))
        n_links = len(data.get("links", []))
        print(
            f"[Python][DataBridge] Export received: {working_fn} "
            f"({n_nodes} nodes, {n_links} links)"
        )

        # Emit so the MainWindow (or any other subscriber) can react
        self.exportDataReceived.emit(data)


# ── Main window ──────────────────────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self, html_path: str, status_note: str = None):
        super().__init__()
        self.html_path = os.path.abspath(html_path)
        self.status_note = status_note
        self.setWindowTitle(f"NMR Viewer — {os.path.basename(html_path)}")
        self.setGeometry(100, 100, 1200, 800)

        # Status bar gives visible feedback when export fires
        self.status = QStatusBar()
        self.setStatusBar(self.status)

        # ── Set up QWebChannel ───────────────────────────────────────────
        self.bridge = DataBridge()
        self.bridge.exportDataReceived.connect(self.on_export_received)

        self.channel = QWebChannel()
        self.channel.registerObject("pyBridge", self.bridge)

        # ── Browser widget ───────────────────────────────────────────────
        self.browser = QWebEngineView()

        # Use the logging page subclass so JS console.log/warn/error calls
        # (including everything in INJECT_JS) actually show up in stdout.
        self.page = LoggingWebEnginePage(self.browser)
        self.browser.setPage(self.page)
        self.page.setWebChannel(self.channel)
        self.page.renderProcessTerminated.connect(self.on_render_process_terminated)

        # ── Auto-grant clipboard permission ───────────────────────────────
        # exportToMnova() calls navigator.clipboard.writeText(). If Qt asks
        # for permission and nothing answers, the request can hang forever
        # (see the writeText timeout patch in INJECT_JS for the safety net
        # regardless — this just makes the real clipboard write succeed
        # too, when the platform supports it).
        if HAS_NEW_PERMISSION_API and hasattr(self.page, "permissionRequested"):
            self.page.permissionRequested.connect(self.on_permission_requested)
        elif hasattr(self.page, "featurePermissionRequested"):
            self.page.featurePermissionRequested.connect(self.on_feature_permission_requested)
        else:
            print("[Python] No permission-request signal found on this PySide6 version — "
                  "relying on the JS-side clipboard timeout guard only.")

        # Newer QtWebEngine blocks file:// pages from fetching remote
        # resources (e.g. the d3.js CDN <script> tag) by default — the old
        # PyQt5 5.15 WebEngine build was more permissive here. Opt back in
        # explicitly so the HTML's external script/stylesheet loads work.
        settings = self.browser.settings()
        settings.setAttribute(
            QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, True
        )
        settings.setAttribute(
            QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls, True
        )

        self.browser.loadFinished.connect(self.on_load_finished)
        self.browser.load(QUrl.fromLocalFile(self.html_path))

        self.setCentralWidget(self.browser)
        self.status.showMessage("Loading…")

    def on_permission_requested(self, permission):
        """Newer PySide6: permission is a QWebEnginePermission object."""
        try:
            permission.grant()
        except Exception as e:
            print(f"[Python] Failed to grant permission: {e}")

    def on_feature_permission_requested(self, origin, feature):
        """Older PySide6: origin (QUrl) + feature (QWebEnginePage.Feature)."""
        self.page.setFeaturePermission(
            origin, feature, QWebEnginePage.PermissionPolicy.PermissionGrantedByUser
        )

    def on_load_finished(self, ok: bool):
        if ok:
            self.status.showMessage(
                self.status_note or "Ready — press Export to send data to Python"
            )
            self.page.runJavaScript(INJECT_JS)
        else:
            self.status.showMessage("ERROR: page failed to load")
            print(f"[Python] ERROR: page failed to load: {self.html_path}")

    def on_render_process_terminated(self, status, exit_code):
        """Called when the Chromium renderer process crashes or is killed.

        Without this, a renderer crash (e.g. the macOS Mach IPC "msg too
        large" failure) leaves the window blank with no Python-visible
        signal at all - the only trace was previously whatever JASON
        happened to capture from stderr.
        """
        status_name = TERMINATION_STATUS_NAMES.get(status, str(status))
        print(
            f"[Python] Renderer process terminated: {status_name} "
            f"(exit code {exit_code}) while loading {self.html_path}"
        )
        self.status.showMessage(
            f"Viewer crashed: {status_name} (exit code {exit_code})"
        )

        log_path = _log_viewer_crash(self.html_path, status_name, exit_code)
        log_note = f"\n\nDetails logged to:\n{log_path}" if log_path else ""

        QMessageBox.critical(
            self,
            "Viewer Crashed",
            "The results viewer's rendering process terminated unexpectedly "
            f"({status_name}, exit code {exit_code}).\n\n"
            "The HTML file itself was saved successfully and can still be "
            f"opened directly in a web browser:\n{self.html_path}{log_note}",
        )

    def on_export_received(self, data: dict):
        """Called on the Python/Qt side when the JavaScript export fires."""
        msg = (
            f"Export received: {data.get('workingFilename', '?')}  |  "
            f"{len(data.get('nodes_now', []))} nodes  |  "
            f"{len(data.get('links', []))} links"
        )
        self.status.showMessage(msg)

        # ── Save to disk (optional) ──────────────────────────────────────
        out_dir = os.path.dirname(self.html_path)
        out_file = os.path.join(out_dir, f"{data.get('workingFilename', 'export')}_py.json")
        with open(out_file, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2)
        print(f"[Python][MainWindow] Saved to {out_file}")


# ── Entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python displayHTML_pyside.py <path_to_html_file>")
        sys.exit(1)

    app = QApplication(sys.argv)
    window = MainWindow(sys.argv[1])
    window.show()
    sys.exit(app.exec())
