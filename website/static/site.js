/* Centrifuger docs — theme switch, mobile nav, copy buttons, TOC scroll-spy
   and the client-side search dialog.  No dependencies. */
(function () {
  'use strict';

  var root = document.documentElement;
  var base = window.CFR_BASE || '';

  /* ------------------------------------------------------------ theme -- */

  var THEMES = ['auto', 'light', 'dark'];
  var themeBtn = document.querySelector('.theme-toggle');
  if (themeBtn) {
    themeBtn.addEventListener('click', function () {
      var current = root.dataset.theme || 'auto';
      var next = THEMES[(THEMES.indexOf(current) + 1) % THEMES.length];
      root.dataset.theme = next;
      themeBtn.title = 'Theme: ' + next + ' (click to change)';
      try { localStorage.setItem('cfr-theme', next); } catch (e) {}
    });
  }

  /* ------------------------------------------------------- mobile nav -- */

  var menuBtn = document.querySelector('.menu-toggle');
  if (menuBtn) {
    menuBtn.addEventListener('click', function () {
      var open = document.body.classList.toggle('nav-open');
      menuBtn.setAttribute('aria-expanded', open ? 'true' : 'false');
    });
    document.addEventListener('click', function (ev) {
      if (!document.body.classList.contains('nav-open')) return;
      if (ev.target.closest('.sidebar') || ev.target.closest('.menu-toggle')) return;
      document.body.classList.remove('nav-open');
      menuBtn.setAttribute('aria-expanded', 'false');
    });
  }

  /* ----------------------------------------------------- copy buttons -- */

  document.querySelectorAll('figure.code [data-copy]').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var code = btn.closest('figure.code').querySelector('code');
      var text = code ? code.textContent : '';
      var done = function () {
        btn.textContent = 'Copied';
        btn.classList.add('done');
        setTimeout(function () {
          btn.textContent = 'Copy';
          btn.classList.remove('done');
        }, 1600);
      };
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(done, fallback);
      } else {
        fallback();
      }
      function fallback() {
        var ta = document.createElement('textarea');
        ta.value = text;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        try { document.execCommand('copy'); done(); } catch (e) {}
        document.body.removeChild(ta);
      }
    });
  });

  /* -------------------------------------------------------- scrollspy -- */

  var tocLinks = Array.prototype.slice.call(document.querySelectorAll('.toc a'));
  if (tocLinks.length) {
    var targets = tocLinks
      .map(function (a) { return document.getElementById(a.hash.slice(1)); })
      .filter(Boolean);

    var sync = function () {
      var top = window.scrollY + 120;
      var current = targets[0];
      targets.forEach(function (el) { if (el.offsetTop <= top) current = el; });
      tocLinks.forEach(function (a) {
        a.classList.toggle('active', current && a.hash === '#' + current.id);
      });
    };
    var ticking = false;
    window.addEventListener('scroll', function () {
      if (ticking) return;
      ticking = true;
      window.requestAnimationFrame(function () { sync(); ticking = false; });
    }, { passive: true });
    sync();
  }

  /* ----------------------------------------------------------- search -- */

  var modal = document.querySelector('.search-modal');
  var input = document.querySelector('.search-input');
  var list = document.querySelector('.search-results');
  var openBtns = document.querySelectorAll('.search-open');
  var index = window.CFR_SEARCH_INDEX || [];
  var selected = 0;

  function openSearch() {
    if (!modal) return;
    modal.hidden = false;
    input.value = '';
    render([]);
    input.focus();
  }

  function closeSearch() {
    if (!modal) return;
    modal.hidden = true;
  }

  openBtns.forEach(function (b) { b.addEventListener('click', openSearch); });

  document.addEventListener('keydown', function (ev) {
    if (!ev.key) return;
    if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === 'k') {
      ev.preventDefault();
      modal && modal.hidden ? openSearch() : closeSearch();
      return;
    }
    if (ev.key === 'Escape') closeSearch();
    if (ev.key === '/' && document.activeElement === document.body) {
      ev.preventDefault();
      openSearch();
    }
  });

  if (modal) {
    modal.addEventListener('click', function (ev) {
      if (ev.target === modal) closeSearch();
    });
  }

  function score(page, terms) {
    var title = page.title.toLowerCase();
    var desc = (page.description || '').toLowerCase();
    var text = (page.text || '').toLowerCase();
    var heads = page.headings.map(function (h) { return h.text.toLowerCase(); }).join(' ');
    var total = 0;
    for (var i = 0; i < terms.length; i++) {
      var t = terms[i];
      var hit = 0;
      if (title.indexOf(t) !== -1) hit += title.indexOf(t) === 0 ? 14 : 10;
      if (heads.indexOf(t) !== -1) hit += 5;
      if (desc.indexOf(t) !== -1) hit += 3;
      if (text.indexOf(t) !== -1) hit += 1;
      if (!hit) return 0;          // every term must appear somewhere
      total += hit;
    }
    return total;
  }

  function snippet(page, term) {
    var text = page.text || page.description || '';
    var at = text.toLowerCase().indexOf(term);
    if (at === -1) return page.description || text.slice(0, 120);
    var start = Math.max(0, at - 45);
    return (start ? '…' : '') + text.slice(start, start + 140).trim() + '…';
  }

  function escapeHtml(s) {
    return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  function render(results) {
    list.innerHTML = '';
    selected = 0;
    if (!results.length) {
      var li = document.createElement('li');
      li.className = 'empty';
      li.textContent = input.value.trim()
        ? 'No results for “' + input.value.trim() + '”.'
        : 'Search page titles, headings and option names.';
      list.appendChild(li);
      return;
    }
    results.forEach(function (r, i) {
      var li = document.createElement('li');
      if (i === 0) li.className = 'sel';
      li.innerHTML =
        '<a href="' + base + r.page.url + '">' +
        '<span class="r-group">' + escapeHtml(r.page.group) + '</span>' +
        '<span class="r-title">' + escapeHtml(r.page.title) + '</span>' +
        '<span class="r-snippet">' + escapeHtml(r.snippet) + '</span></a>';
      list.appendChild(li);
    });
  }

  function search(query) {
    var terms = query.toLowerCase().split(/\s+/).filter(Boolean);
    if (!terms.length) return [];
    return index
      .map(function (page) { return { page: page, s: score(page, terms) }; })
      .filter(function (r) { return r.s > 0; })
      .sort(function (a, b) { return b.s - a.s; })
      .slice(0, 8)
      .map(function (r) {
        r.snippet = snippet(r.page, terms[0]);
        return r;
      });
  }

  if (input) {
    input.addEventListener('input', function () { render(search(input.value)); });
    input.addEventListener('keydown', function (ev) {
      var items = list.querySelectorAll('li a');
      if (!items.length) return;
      if (ev.key === 'ArrowDown' || ev.key === 'ArrowUp') {
        ev.preventDefault();
        list.children[selected].classList.remove('sel');
        selected = (selected + (ev.key === 'ArrowDown' ? 1 : items.length - 1)) % items.length;
        list.children[selected].classList.add('sel');
        list.children[selected].scrollIntoView({ block: 'nearest' });
      } else if (ev.key === 'Enter') {
        ev.preventDefault();
        items[selected].click();
      }
    });
  }
})();
