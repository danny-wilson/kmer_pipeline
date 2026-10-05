"""report_assets.py: self-contained HTML reports (D12). The style sheet and script are written
into each report and the figures embedded as data URIs, so a report can be opened on its own:
moved, attached to an e-mail, or shown inside JupyterLab, whose viewer cannot fetch files next
to the report. The content is visible without the script (JupyterLab runs it only once the
report is trusted): the script only turns the figures into slide shows."""
import base64
import os


def head_style(src_dir):
    with open(os.path.join(src_dir, "report.css")) as fh:
        return "  <style>\n" + fh.read() + "  </style>"


def foot_script(src_dir):
    with open(os.path.join(src_dir, "report.js")) as fh:
        return "<script>\n" + fh.read() + "</script>"


def img(path, attrs):
    """An <img> tag with the figure embedded, or "" if the figure was not drawn."""
    if not path or not os.path.isfile(path):
        return ""
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return '<img src="data:image/png;base64,' + data + '" ' + attrs + '>'
