"""Self-contained HTML reports (D12)."""
import base64
import os
import re

from conftest import REPO_DIR, SCRIPTS_DIR

import report_assets

ASSETS = next(d for d in (SCRIPTS_DIR, REPO_DIR) if os.path.exists(os.path.join(d, "report.css")))
PNG = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg==")


def test_img(tmp_path):
    f = tmp_path / "fig.png"
    f.write_bytes(PNG)
    tag = report_assets.img(str(f), 'class="center"')
    m = re.fullmatch(r'<img src="data:image/png;base64,([^"]+)" class="center">', tag)
    assert m and base64.b64decode(m.group(1)) == PNG
    assert report_assets.img(str(tmp_path / "missing.png"), 'class="x"') == ""
    assert report_assets.img(None, 'class="x"') == ""


def test_style_and_script_inline():
    assert report_assets.head_style(ASSETS).startswith("  <style>") and ".mySlides" in report_assets.head_style(ASSETS)
    assert "function showSlides" in report_assets.foot_script(ASSETS)


def test_figures_visible_without_the_script():
    css = open(os.path.join(ASSETS, "report.css")).read()
    for cls in ("mySlides", "toggled", "untoggled"):
        block = re.search(r"\." + cls + r" \{([^}]*)\}", css).group(1)
        assert "display: block" in block, cls


def check_report(path):
    """Problems with a generated report: any src that is not embedded, any link that is not to
    an existing report next to it."""
    html = open(path).read()
    problems = [s for s in re.findall(r'src="([^"]*)"', html) if not s.startswith("data:image/png;base64,")]
    problems += [s for s in re.findall(r"src='([^']*)'", html)]
    folder = os.path.dirname(path)
    for h in re.findall(r"href='([^']*)'", html) + re.findall(r'href="([^"]*)"', html):
        if not (h.endswith(".html") and "/" not in h and os.path.exists(os.path.join(folder, h))):
            problems.append(h)
    return problems
