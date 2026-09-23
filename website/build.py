#!/usr/bin/env python3
"""Static site generator for the Centrifuger documentation website.

No third-party dependencies: it only needs a Python 3.6+ interpreter.  Markdown
pages under ``website/content`` are rendered into a self-contained static site
under ``docs/`` that GitHub Pages can serve directly.

    python3 website/build.py              # build into ./docs
    python3 website/build.py -o /tmp/out  # build somewhere else

Links between pages are written site-root-absolute in the Markdown sources
(``/guides/classification/``) and are rewritten to *relative* URLs at build
time, so the generated site works unchanged at https://user.github.io/repo/
and at a domain root.  It must be served over HTTP rather than opened from
disk, because the clean URLs resolve to directories.
"""

from __future__ import print_function

import argparse
import json
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONTENT = os.path.join(HERE, "content")
STATIC = os.path.join(HERE, "static")

SITE_TITLE = "Centrifuger"
SITE_TAGLINE = "Taxonomic classification on a compressed genome database"
REPO_URL = "https://github.com/mourisl/centrifuger"
EDIT_BASE = "https://github.com/mourisl/centrifuger/edit/master/website/content/"

# --------------------------------------------------------------------------
# Site navigation.  Each entry is (group label, [(page label, slug), ...]).
# The slug is the path of the Markdown file relative to content/, without the
# .md extension, and also the URL of the generated page.
# --------------------------------------------------------------------------
NAV = [
    ("Getting started", [
        ("Introduction", "getting-started/introduction"),
        ("Installation", "getting-started/installation"),
        ("Quick start", "getting-started/quick-start"),
    ]),
    ("Guides", [
        ("Building an index", "guides/building-an-index"),
        ("Pre-built indexes", "guides/prebuilt-indexes"),
        ("Classifying reads", "guides/classification"),
        ("Abundance quantification", "guides/quantification"),
        ("Single-cell and barcoded data", "guides/single-cell"),
        ("Protein classification", "guides/protein"),
    ]),
    ("Reference", [
        ("Command-line interface", "reference/cli"),
        ("Output formats", "reference/output-formats"),
        ("FAQ and troubleshooting", "reference/faq"),
    ]),
    ("About", [
        ("Citation and support", "about/citation"),
    ]),
]

# Flat, ordered list of slugs used for the previous/next footer links.
FLAT_NAV = [(label, slug) for _, items in NAV for label, slug in items]


# ==========================================================================
# Markdown -> HTML
# ==========================================================================

_ESCAPES = (("&", "&amp;"), ("<", "&lt;"), (">", "&gt;"))


def esc(text):
    for a, b in _ESCAPES:
        text = text.replace(a, b)
    return text


def slugify(text):
    text = re.sub(r"<[^>]+>", "", text)
    text = text.replace("&amp;", "and").replace("&lt;", "").replace("&gt;", "")
    text = text.lower()
    text = re.sub(r"[^a-z0-9]+", "-", text)
    return text.strip("-") or "section"


class Renderer(object):
    """Renders the Markdown subset used by this site.

    Supported: ATX headings, fenced code blocks, paragraphs, nested bullet and
    ordered lists, GFM tables, blockquotes, ``:::type`` admonitions, thematic
    breaks, raw HTML blocks, and the inline constructs handled by ``inline()``.
    """

    def __init__(self, prefix):
        # ``prefix`` turns a site-root-absolute URL into a relative one.
        self.prefix = prefix
        self.headings = []
        self.seen_ids = {}

    # ---------------- inline ----------------

    def _unique_id(self, base):
        if base not in self.seen_ids:
            self.seen_ids[base] = 1
            return base
        self.seen_ids[base] += 1
        return "%s-%d" % (base, self.seen_ids[base])

    def _url(self, url):
        if url.startswith("/"):
            return self.prefix + url.lstrip("/")
        return url

    def _link(self, match):
        text, url = match.group(1), match.group(2)
        extra = ""
        if url.startswith("http"):
            extra = ' target="_blank" rel="noopener"'
        return '<a href="%s"%s>%s</a>' % (self._url(url), extra, text)

    def inline(self, text):
        """Inline formatting.

        Code spans are pulled out into placeholders *before* anything else
        runs, and restored at the very end.  Doing it that way keeps their
        contents safe from the emphasis and link rules while still allowing
        constructs that wrap a code span, such as **`numReads`** or
        [`centrifuger-quant`](/guides/quantification/).
        """
        spans = []

        def stash(match):
            spans.append(match.group(1))
            return "\x00%d\x00" % (len(spans) - 1)

        text = re.sub(r"`([^`]+)`", stash, text)

        text = esc(text)
        text = re.sub(r"!\[([^\]]*)\]\(([^)\s]+)\)",
                      lambda m: '<img src="%s" alt="%s">' % (self._url(m.group(2)), m.group(1)),
                      text)
        text = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", self._link, text)
        text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
        text = re.sub(r"(?<!\*)\*([^*\s][^*]*?)\*(?!\*)", r"<em>\1</em>", text)

        return re.sub(r"\x00(\d+)\x00",
                      lambda m: "<code>%s</code>" % esc(spans[int(m.group(1))]),
                      text)

    def _rewrite_raw_urls(self, html):
        """Raw HTML blocks bypass inline(), so site-root-absolute href/src
        attributes inside them still need turning into relative URLs."""
        return re.sub(r'(href|src)="/([^"]*)"',
                      lambda m: '%s="%s%s"' % (m.group(1), self.prefix, m.group(2)),
                      html)

    # ---------------- blocks ----------------

    def render(self, lines):
        html = []
        i = 0
        n = len(lines)
        while i < n:
            line = lines[i]
            stripped = line.strip()

            if not stripped:
                i += 1
                continue

            # Fenced code block
            if stripped.startswith("```"):
                lang = stripped[3:].strip()
                body = []
                i += 1
                while i < n and not lines[i].strip().startswith("```"):
                    body.append(lines[i])
                    i += 1
                i += 1
                html.append(self._code_block(body, lang))
                continue

            # Admonition
            if stripped.startswith(":::"):
                header = stripped[3:].strip()
                kind, _, title = header.partition(" ")
                body = []
                i += 1
                while i < n and not lines[i].strip() == ":::":
                    body.append(lines[i])
                    i += 1
                i += 1
                html.append(self._admonition(kind or "note", title.strip(), body))
                continue

            # Thematic break
            if re.match(r"^(\*\s*){3,}$|^(-\s*){3,}$", stripped):
                html.append("<hr>")
                i += 1
                continue

            # Heading
            m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
            if m:
                level = len(m.group(1))
                text = self.inline(m.group(2).strip())
                hid = self._unique_id(slugify(m.group(2)))
                if level in (2, 3):
                    self.headings.append((level, re.sub(r"<[^>]+>", "", text), hid))
                html.append(
                    '<h%d id="%s">%s<a class="anchor" href="#%s" aria-label="Link to this section">#</a></h%d>'
                    % (level, hid, text, hid, level))
                i += 1
                continue

            # Raw HTML block: starts at column 0 with a tag, ends at a blank line
            if line.startswith("<"):
                body = []
                while i < n and lines[i].strip():
                    body.append(lines[i])
                    i += 1
                html.append(self._rewrite_raw_urls("\n".join(body)))
                continue

            # Table
            if "|" in stripped and i + 1 < n and re.match(r"^\s*\|?[\s:|-]+\|[\s:|-]*$", lines[i + 1]):
                block = []
                while i < n and "|" in lines[i] and lines[i].strip():
                    block.append(lines[i])
                    i += 1
                html.append(self._table(block))
                continue

            # Blockquote
            if stripped.startswith(">"):
                body = []
                while i < n and lines[i].strip().startswith(">"):
                    body.append(re.sub(r"^\s*>\s?", "", lines[i]))
                    i += 1
                html.append("<blockquote>%s</blockquote>" % self.render(body))
                continue

            # List
            if re.match(r"^\s*([-*]|\d+\.)\s+", line):
                block = []
                while i < n and (lines[i].strip() and
                                 (re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]) or lines[i].startswith(("  ", "\t")))):
                    block.append(lines[i])
                    i += 1
                html.append(self._list(block))
                continue

            # Paragraph
            para = []
            while i < n and lines[i].strip() and not lines[i].strip().startswith(("```", ":::", "#", ">")) \
                    and not re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]) and not lines[i].startswith("<"):
                para.append(lines[i].strip())
                i += 1
            if para:
                html.append("<p>%s</p>" % self.inline(" ".join(para)))
            else:
                i += 1

        return "\n".join(html)

    def _code_block(self, body, lang):
        # Trim a uniform leading indent so copy/paste gives clean commands.
        while body and not body[0].strip():
            body.pop(0)
        while body and not body[-1].strip():
            body.pop()
        filled = [l for l in body if l.strip()]
        if filled:
            dedent = min(len(l) - len(l.lstrip(" ")) for l in filled)
            body = [l[dedent:] if l.strip() else "" for l in body]
        code = esc("\n".join(body))
        label = lang or "text"
        return ('<figure class="code">'
                '<figcaption><span class="lang">%s</span>'
                '<button class="copy" type="button" data-copy aria-label="Copy code">Copy</button>'
                '</figcaption>'
                '<pre><code class="language-%s">%s</code></pre></figure>'
                % (esc(label), esc(label), code))

    def _admonition(self, kind, title, body):
        kind = kind.lower()
        if kind not in ("note", "tip", "caution", "danger"):
            kind = "note"
        if not title:
            title = {"note": "Note", "tip": "Tip",
                     "caution": "Caution", "danger": "Danger"}[kind]
        return ('<aside class="admonition %s"><p class="admonition-title">%s</p>%s</aside>'
                % (kind, esc(title), self.render(body)))

    def _table(self, block):
        rows = []
        for raw in block:
            cells = raw.strip().strip("|").split("|")
            rows.append([c.strip() for c in cells])
        header, body = rows[0], rows[2:]
        out = ['<div class="table-wrap"><table><thead><tr>']
        for cell in header:
            out.append("<th>%s</th>" % self.inline(cell))
        out.append("</tr></thead><tbody>")
        for row in body:
            # A row whose only non-empty cell is bold acts as a section divider.
            filled = [c for c in row if c]
            if len(filled) == 1 and filled[0].startswith("**"):
                out.append('<tr class="group"><td colspan="%d">%s</td></tr>'
                           % (len(header), self.inline(filled[0])))
                continue
            out.append("<tr>")
            for cell in row[:len(header)]:
                out.append("<td>%s</td>" % self.inline(cell))
            out.append("</tr>")
        out.append("</tbody></table></div>")
        return "".join(out)

    def _list(self, block):
        def indent_of(s):
            return len(s) - len(s.lstrip(" "))

        def build(items, base):
            ordered = bool(re.match(r"^\s*\d+\.\s+", items[0]))
            out = ["<ol>" if ordered else "<ul>"]
            idx = 0
            while idx < len(items):
                line = items[idx]
                m = re.match(r"^\s*(?:[-*]|\d+\.)\s+(.*)$", line)
                if m is None:
                    idx += 1
                    continue
                text = m.group(1)
                child = []
                idx += 1
                while idx < len(items) and indent_of(items[idx]) > base:
                    child.append(items[idx])
                    idx += 1
                inner = self.inline(text)
                if child:
                    dedent = min(indent_of(c) for c in child)
                    inner += build([c[dedent:] for c in child], 0)
                out.append("<li>%s</li>" % inner)
            out.append("</ol>" if ordered else "</ul>")
            return "".join(out)

        base = min(indent_of(l) for l in block if l.strip())
        return build([l[base:] for l in block], 0)


# ==========================================================================
# Page loading
# ==========================================================================

def read_page(slug):
    path = os.path.join(CONTENT, slug + ".md")
    with open(path, "r") as fh:
        raw = fh.read()
    meta = {}
    if raw.startswith("---"):
        _, front, raw = raw.split("---", 2)
        for line in front.strip().splitlines():
            key, _, value = line.partition(":")
            meta[key.strip()] = value.strip()
    return meta, raw.strip("\n").splitlines()


def prefix_for(slug):
    """Relative path back to the site root from the page's directory."""
    if slug == "index":
        return ""
    return "../" * (slug.count("/") + 1)


def url_for(slug):
    return "index.html" if slug == "index" else slug + "/"


# ==========================================================================
# Templates
# ==========================================================================

GITHUB_ICON = (
    '<svg viewBox="0 0 16 16" aria-hidden="true" width="18" height="18">'
    '<path fill="currentColor" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 '
    '0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 '
    '1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 '
    '0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82a7.4 7.4 0 0 1 2-.27c.68 0 1.36.09 2 .27 '
    '1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 '
    '3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8Z"/>'
    '</svg>'
)


def render_sidebar(prefix, active_slug):
    out = ['<nav class="sidebar-nav" aria-label="Documentation">']
    for group, items in NAV:
        out.append('<p class="nav-group">%s</p><ul>' % esc(group))
        for label, slug in items:
            cls = ' class="active"' if slug == active_slug else ""
            aria = ' aria-current="page"' if slug == active_slug else ""
            out.append('<li><a href="%s%s"%s%s>%s</a></li>'
                       % (prefix, url_for(slug), cls, aria, esc(label)))
        out.append("</ul>")
    out.append("</nav>")
    return "".join(out)


def render_toc(headings):
    if len(headings) < 2:
        return ""
    out = ['<nav class="toc" aria-label="On this page"><p class="toc-title">On this page</p><ul>']
    for level, text, hid in headings:
        out.append('<li class="lvl%d"><a href="#%s">%s</a></li>' % (level, hid, text))
    out.append("</ul></nav>")
    return "".join(out)


def render_pager(prefix, slug):
    slugs = [s for _, s in FLAT_NAV]
    if slug not in slugs:
        return ""
    idx = slugs.index(slug)
    out = ['<nav class="pager" aria-label="Page navigation">']
    if idx > 0:
        label, target = FLAT_NAV[idx - 1]
        out.append('<a class="prev" href="%s%s"><span>Previous</span>%s</a>'
                   % (prefix, url_for(target), esc(label)))
    else:
        out.append("<span></span>")
    if idx < len(FLAT_NAV) - 1:
        label, target = FLAT_NAV[idx + 1]
        out.append('<a class="next" href="%s%s"><span>Next</span>%s</a>'
                   % (prefix, url_for(target), esc(label)))
    else:
        out.append("<span></span>")
    out.append("</nav>")
    return "".join(out)


PAGE = u"""<!DOCTYPE html>
<html lang="en" data-theme="auto">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<meta name="description" content="{description}">
<meta property="og:title" content="{title}">
<meta property="og:description" content="{description}">
<meta property="og:type" content="website">
<link rel="stylesheet" href="{prefix}assets/style.css">
<script>
// Applied before first paint so the page never flashes the wrong theme.
(function () {{
  try {{
    var t = localStorage.getItem('cfr-theme') || 'auto';
    document.documentElement.dataset.theme = t;
  }} catch (e) {{}}
}})();
</script>
</head>
<body class="{body_class}">
<a class="skip" href="#main">Skip to content</a>

<header class="site-header">
  <button class="menu-toggle" type="button" aria-label="Toggle navigation" aria-expanded="false">
    <span></span><span></span><span></span>
  </button>
  <a class="brand" href="{prefix}index.html">Centrifuger</a>
  <div class="header-spacer"></div>
  <button class="search-open" type="button" aria-label="Search">
    <svg viewBox="0 0 20 20" width="16" height="16" aria-hidden="true"><path fill="none" stroke="currentColor" stroke-width="2" d="M13 13l4 4M8.5 14a5.5 5.5 0 1 1 0-11 5.5 5.5 0 0 1 0 11Z"/></svg>
    <span class="search-label">Search</span><kbd>Ctrl K</kbd>
  </button>
  <button class="theme-toggle" type="button" aria-label="Switch theme" title="Switch between light, dark and system theme">
    <svg class="icon-sun" viewBox="0 0 24 24" width="18" height="18" aria-hidden="true"><circle cx="12" cy="12" r="4.5" fill="currentColor"/><g stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M12 2v2.5M12 19.5V22M2 12h2.5M19.5 12H22M4.9 4.9l1.8 1.8M17.3 17.3l1.8 1.8M19.1 4.9l-1.8 1.8M6.7 17.3l-1.8 1.8"/></g></svg>
    <svg class="icon-moon" viewBox="0 0 24 24" width="18" height="18" aria-hidden="true"><path fill="currentColor" d="M20 14.5A8.5 8.5 0 0 1 9.5 4a8.5 8.5 0 1 0 10.5 10.5Z"/></svg>
  </button>
  <a class="gh" href="{repo}" target="_blank" rel="noopener" aria-label="Centrifuger on GitHub">{github}</a>
</header>

<div class="layout">
  <aside class="sidebar">{sidebar}</aside>
  <div class="content-col">
    <main id="main">
{article}
    </main>
    <footer class="site-footer">
      <p>Centrifuger is developed by <a href="https://github.com/mourisl" target="_blank" rel="noopener">Li Song</a>
         and colleagues. Copyright &copy; 2023&ndash;present, Li Song.</p>
      <p><a href="{repo}" target="_blank" rel="noopener">Source code</a> &middot;
         <a href="{repo}/issues" target="_blank" rel="noopener">Report an issue</a> &middot;
         <a href="https://doi.org/10.1186/s13059-024-03244-4" target="_blank" rel="noopener">Genome Biology paper</a></p>
    </footer>
  </div>
  {aside}
</div>

<div class="search-modal" hidden>
  <div class="search-panel" role="dialog" aria-modal="true" aria-label="Search the documentation">
    <input type="search" class="search-input" placeholder="Search the documentation&hellip;" autocomplete="off" spellcheck="false">
    <ul class="search-results"></ul>
    <p class="search-hint">Enter to open &middot; Esc to close</p>
  </div>
</div>

<script>window.CFR_BASE = "{prefix}";</script>
<script src="{prefix}assets/search.js" defer></script>
<script src="{prefix}assets/site.js" defer></script>
</body>
</html>
"""


def build_page(slug, out_root):
    meta, lines = read_page(slug)
    prefix = prefix_for(slug)
    renderer = Renderer(prefix)
    body = renderer.render(lines)

    is_home = slug == "index"
    title = meta.get("title", SITE_TITLE)
    description = meta.get("description", SITE_TAGLINE)

    if is_home:
        article = body
        aside = ""
    else:
        edit = ('<a class="edit-link" href="%s%s.md" target="_blank" rel="noopener">Edit this page on GitHub</a>'
                % (EDIT_BASE, slug))
        article = ('<article class="prose">\n<div class="page-head"><p class="eyebrow">%s</p>'
                   '<h1>%s</h1><p class="lede">%s</p></div>\n%s\n%s\n%s\n</article>'
                   % (esc(group_of(slug)), esc(title), esc(description), body,
                      render_pager(prefix, slug), edit))
        aside = '<aside class="rightbar">%s</aside>' % render_toc(renderer.headings)

    html = PAGE.format(
        title=esc(title if is_home else "%s | %s" % (title, SITE_TITLE)),
        description=esc(description),
        prefix=prefix,
        body_class="home" if is_home else "doc",
        github=GITHUB_ICON,
        repo=REPO_URL,
        sidebar=render_sidebar(prefix, slug),
        article=article,
        aside=aside,
    )

    if is_home:
        dest = os.path.join(out_root, "index.html")
    else:
        dest = os.path.join(out_root, slug, "index.html")
    ensure_dir(os.path.dirname(dest))
    with open(dest, "w") as fh:
        fh.write(html)

    return {
        "title": title,
        "description": description,
        "url": url_for(slug),
        "group": group_of(slug),
        "headings": [{"text": h[1], "id": h[2]} for h in renderer.headings],
        "text": plain_text(lines),
    }


def group_of(slug):
    for group, items in NAV:
        for _, s in items:
            if s == slug:
                return group
    return SITE_TITLE


def plain_text(lines):
    """Flatten a page to searchable plain text (code blocks included, since
    option names such as --build-mem are exactly what people search for)."""
    out = []
    for line in lines:
        line = line.strip()
        if line.startswith(":::") or line.startswith("```") or line.startswith("<"):
            continue
        line = re.sub(r"[#>*`|]", " ", line)
        line = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", line)
        out.append(line)
    text = re.sub(r"\s+", " ", " ".join(out))
    return text[:12000]


def ensure_dir(path):
    if path and not os.path.isdir(path):
        os.makedirs(path)


def main():
    parser = argparse.ArgumentParser(description="Build the Centrifuger documentation site.")
    parser.add_argument("-o", "--output", default=os.path.join(os.path.dirname(HERE), "docs"),
                        help="output directory (default: ../docs)")
    parser.add_argument("--site-base", default="/centrifuger/",
                        help="absolute URL path the site is served from; used only by "
                             "404.html, which GitHub Pages serves from arbitrary URLs "
                             "and so cannot use relative links (default: /centrifuger/)")
    args = parser.parse_args()
    out_root = os.path.abspath(args.output)

    if os.path.isdir(out_root):
        shutil.rmtree(out_root)
    ensure_dir(out_root)

    slugs = ["index"] + [s for _, s in FLAT_NAV]
    index = []
    for slug in slugs:
        record = build_page(slug, out_root)
        if slug != "index":
            index.append(record)
        print("  rendered %s" % url_for(slug))

    assets = os.path.join(out_root, "assets")
    ensure_dir(assets)
    for name in sorted(os.listdir(STATIC)):
        shutil.copyfile(os.path.join(STATIC, name), os.path.join(assets, name))
        print("  copied  assets/%s" % name)

    with open(os.path.join(assets, "search.js"), "w") as fh:
        fh.write("window.CFR_SEARCH_INDEX = %s;\n"
                 % json.dumps(index, separators=(",", ":")))

    # Tell GitHub Pages not to run the output through Jekyll.
    open(os.path.join(out_root, ".nojekyll"), "w").close()

    # A 404 page that keeps the site chrome.
    with open(os.path.join(out_root, "404.html"), "w") as fh:
        fh.write(PAGE.format(
            title="Page not found | Centrifuger",
            description="The requested page does not exist.",
            prefix=args.site_base,
            body_class="doc",
            github=GITHUB_ICON, repo=REPO_URL,
            sidebar=render_sidebar(args.site_base, ""),
            article=('<article class="prose"><div class="page-head"><p class="eyebrow">Error</p>'
                     '<h1>Page not found</h1><p class="lede">That page does not exist. '
                     'Try the navigation on the left, or search with Ctrl&nbsp;K.</p></div></article>'),
            aside=""))

    print("\nSite written to %s" % out_root)
    print("Preview with:  (cd %s && python3 -m http.server 8000)" % out_root)


if __name__ == "__main__":
    sys.exit(main())
