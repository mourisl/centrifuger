# Centrifuger documentation website

The source for the Centrifuger documentation site. Markdown pages in
`content/` are rendered by `build.py` into a static site in `../docs/`,
which GitHub Pages serves as-is.

```
website/
├── build.py            generator (Python 3.6+, standard library only)
├── content/            the pages, as Markdown
│   ├── index.md            landing page
│   ├── getting-started/
│   ├── guides/
│   ├── reference/
│   └── about/
└── static/             style.css, site.js
```

## Building

```bash
python3 website/build.py
```

That rewrites `docs/` from scratch. There is nothing to install: the
generator uses only the Python standard library, and the site it produces
loads no external scripts, fonts or stylesheets.

Preview it locally:

```bash
(cd docs && python3 -m http.server 8000)
# then open http://localhost:8000/
```

Every internal link is relative, so the site works at any URL prefix
without configuration. It does need to be served over HTTP, though:
links point at clean URLs such as `guides/classification/`, and a browser
opening those over `file://` shows a directory listing rather than the
page inside it. The one-line server above is enough.

## Deploying to GitHub Pages

Commit the generated `docs/` directory, then in the repository go to
**Settings → Pages** and set the source to **Deploy from a branch**, with
branch `master` (or `main`) and folder `/docs`. The site appears at
`https://<user>.github.io/<repo>/`.

`docs/.nojekyll` is generated automatically so that Pages serves the files
untouched instead of running them through Jekyll.

Alternatively, `.github/workflows/docs.yml` builds and publishes the site
with GitHub Actions on every push, which means contributors only have to
edit Markdown. Enable it by setting **Settings → Pages → Source** to
**GitHub Actions**. Use one mechanism or the other, not both.

### A different URL prefix

Page-to-page links are relative and need no configuration. The only file
that must know the absolute URL is `404.html`, which GitHub Pages serves
from arbitrary paths:

```bash
python3 website/build.py --site-base /centrifuger/
```

## Editing content

Each page is a Markdown file with a short front matter block:

```markdown
---
title: Classifying reads
description: One sentence, shown under the page title and in search results.
---
```

The page's position in the sidebar, and the previous/next links at the
bottom, come from the `NAV` list near the top of `build.py`. Adding a page
means creating the Markdown file and adding one line to `NAV`.

### Supported Markdown

Headings, paragraphs, bullet and numbered lists, GitHub-flavoured tables,
blockquotes, fenced code blocks, horizontal rules, and the usual inline
formatting (`code`, **bold**, *italic*, links, images). Raw HTML blocks pass
through untouched, which is how the landing page builds its hero and card
grid; such a block must start at column 0 and contain no blank lines.

Two conventions are specific to this site:

**Links between pages** are written site-root-absolute and rewritten to
relative URLs at build time:

```markdown
See [building an index](/guides/building-an-index/).
```

**Admonitions** use a fenced syntax, with an optional title:

```markdown
:::caution The quotes are required
Without them the shell expands the pattern before Centrifuger sees it.
:::
```

Available types: `note`, `tip`, `caution`, `danger`.

### Things the generator does not do

Fenced code blocks and tables must be written at the top level, not nested
inside a list item. Reference-style links (`[text][ref]`) and footnotes are
not supported. These limits keep the generator small and dependency-free;
if a page needs more, it can drop into raw HTML.

## Branding

Centrifuger has no official logo, so the site does not use one: the header
is the project name as text, and there is no favicon. Please do not add an
invented mark.

## What the site includes

- Light, dark and system themes, remembered per visitor, applied before
  first paint so there is no flash of the wrong theme.
- Client-side search over page titles, headings and body text, opened with
  <kbd>Ctrl</kbd>+<kbd>K</kbd> or <kbd>/</kbd>. The index is generated into
  `docs/assets/search.js` at build time.
- An on-page table of contents with scroll tracking, copy buttons on code
  blocks, previous/next navigation, and an off-canvas sidebar on narrow
  screens.
