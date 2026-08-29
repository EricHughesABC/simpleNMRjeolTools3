"""
nmr_viewer.py
─────────────
Loads a simpleNMR D3 HTML file in a QWebEngineView and intercepts the
"Export" button so the JSON data is sent back to Python via QWebChannel
instead of being saved to disk.

No modifications to the HTML file are required.

Usage:
    python nmr_viewer.py path/to/Cytochalasin-B-kate_d3.html

Requirements:
    pip install PyQt5 PyQtWebEngine
"""

import sys
import os
import json

# ── QWebEngineWidgets MUST be imported before QApplication ──────────────────
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWebChannel import QWebChannel
from PyQt5.QtWidgets import QApplication, QMainWindow, QStatusBar
from PyQt5.QtCore import QObject, QUrl, pyqtSignal, pyqtSlot

# Helps prevent sandbox errors on macOS
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--no-sandbox"

# ── JavaScript injected after the page has finished loading ─────────────────
# Strategy:
#   1. Dynamically load qwebchannel.js from Qt's built-in resource path.
#   2. Once loaded, open the channel and store window.pyBridge.
#   3. Patch HTMLAnchorElement.prototype.click so that any programmatic
#      anchor.click() on a data:text/json URL is intercepted and the JSON
#      is sent to Python instead of triggering a browser download.
#      If pyBridge is not available (standalone browser use), the original
#      click behaviour is preserved as a fallback.
INJECT_JS = """
(function () {
    // ── 1. Load qwebchannel.js from Qt's built-in resource ──────────────
    var script = document.createElement('script');
    script.src = 'qrc:///qtwebchannel/qwebchannel.js';

    script.onload = function () {

        // ── 2. Open the channel Qt registered as "pyBridge" ─────────────
        new QWebChannel(qt.webChannelTransport, function (channel) {
            window.pyBridge = channel.objects.pyBridge;
            console.log('[nmr_viewer] QWebChannel ready, pyBridge connected.');
        });
    };

    document.head.appendChild(script);

    // ── 3. Patch anchor clicks to intercept JSON data-URI downloads ──────
    //    We patch the prototype once, before any user interaction occurs.
    //    exportToMnova() calls anchor.click() programmatically, so this
    //    fires reliably.
    var _origClick = HTMLAnchorElement.prototype.click;

    HTMLAnchorElement.prototype.click = function () {
        var href = this.href || '';

        if (href.startsWith('data:text/json') && window.pyBridge) {
            // Decode the JSON from the data URI
            var encoded = href.replace(/^data:text\/json;charset=utf-8,/, '');
            var jsonStr = decodeURIComponent(encoded);

            // Send to Python — this calls DataBridge.receiveExportData()
            window.pyBridge.receiveExportData(jsonStr);
            console.log('[nmr_viewer] Export intercepted and sent to Python.');

            // Do NOT call _origClick — skip the file download
            return;
        }

        // Anything else (non-JSON links, etc.) behaves normally
        _origClick.call(this);
    };

    console.log('[nmr_viewer] Anchor-click patch applied.');
})();
"""


# ── Python-side bridge object exposed to JavaScript ─────────────────────────
class DataBridge(QObject):
    """
    Registered with QWebChannel under the name 'pyBridge'.
    JavaScript calls  window.pyBridge.receiveExportData(jsonString)
    which triggers the receiveExportData slot below.
    """

    # Emitted with the parsed dict whenever export data arrives
    exportDataReceived = pyqtSignal(dict)

    @pyqtSlot(str)
    def receiveExportData(self, json_str: str):
        """Called from JavaScript when the Export button is pressed."""
        try:
            data = json.loads(json_str)
        except json.JSONDecodeError as e:
            print(f"[DataBridge] JSON decode error: {e}")
            return

        print(f"[DataBridge] Export received!  Keys: {list(data.keys())}")

        # ── Do whatever you need with the data here ──────────────────────
        # Examples:
        #   • Save to a specific location:
        #       out_path = Path(data.get('workingDirectory', '.')) / 'export.json'
        #       out_path.write_text(json_str)
        #   • Pass to another part of your application:
        #       self.exportDataReceived.emit(data)
        #   • Print a summary:
        working_fn  = data.get('workingFilename', 'unknown')
        n_nodes     = len(data.get('nodes_now', []))
        n_links     = len(data.get('links', []))
        print(f"           Molecule : {working_fn}")
        print(f"           Nodes    : {n_nodes}")
        print(f"           Links    : {n_links}")

        # Emit so the MainWindow (or any other subscriber) can react
        self.exportDataReceived.emit(data)


# ── Main window ──────────────────────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self, html_path: str):
        super().__init__()
        self.html_path = os.path.abspath(html_path)
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
        self.browser.page().setWebChannel(self.channel)
        self.browser.loadFinished.connect(self.on_load_finished)
        self.browser.load(QUrl.fromLocalFile(self.html_path))

        self.setCentralWidget(self.browser)
        self.status.showMessage("Loading…")

    def on_load_finished(self, ok: bool):
        if ok:
            self.status.showMessage("Ready — press Export to send data to Python")
            # Inject QWebChannel setup + anchor-click patch
            self.browser.page().runJavaScript(INJECT_JS)
        else:
            self.status.showMessage("ERROR: page failed to load")

    def on_export_received(self, data: dict):
        """Called on the Python/Qt side when the JavaScript export fires."""
        msg = (
            f"Export received: {data.get('workingFilename', '?')}  |  "
            f"{len(data.get('nodes_now', []))} nodes  |  "
            f"{len(data.get('links', []))} links"
        )
        self.status.showMessage(msg)

        # ── Save to disk (optional) ──────────────────────────────────────
        out_dir  = os.path.dirname(self.html_path)
        out_file = os.path.join(out_dir, f"{data.get('workingFilename', 'export')}_py.json")
        with open(out_file, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2)
        print(f"[MainWindow] Saved to {out_file}")


# ── Entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python nmr_viewer.py <path_to_html_file>")
        sys.exit(1)

    app = QApplication(sys.argv)
    window = MainWindow(sys.argv[1])
    window.show()
    sys.exit(app.exec_())