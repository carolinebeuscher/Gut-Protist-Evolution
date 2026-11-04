#!/usr/bin/env python3
"""
Fetches the KEGG pathway hierarchy from https://www.genome.jp/kegg/pathway.html
and writes a TSV with four columns:

    Pathway   Group                           Supergroup   Annot
    00010     Carbohydrate metabolism         Metabolism   Glycolysis / Gluconeogenesis

Usage:
    python3 build_kegg_pathway_hierarchy.py --output kegg_pathway_hierarchy.tsv

Install dependencies:
    pip install requests beautifulsoup4
"""

import argparse
import re
import sys

try:
    import requests
    from bs4 import BeautifulSoup
except ImportError:
    sys.exit("ERROR: run: pip install requests beautifulsoup4")

URL = "https://www.genome.jp/kegg/pathway.html"


def fetch_page(url):
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text


def parse_hierarchy(html):
    """
    Parse the KEGG pathway page into rows of
    (pathway_id, group, supergroup, annot).

    Page structure:
      <h4>  → top-level section number + name  → Supergroup
              e.g. "1. Metabolism", "2. Genetic Information Processing"
      <b>   → subsection number + name          → Group
              e.g. "1.1 Carbohydrate metabolism"
      text  → pathway number + annotation        → Pathway + Annot
              e.g. "00010 M N\nGlycolysis / Gluconeogenesis"
    """
    soup = BeautifulSoup(html, "html.parser")
    rows = []

    # We walk the content div linearly, tracking current supergroup and group
    supergroup = ""
    group      = ""

    # Find the main content area — everything after "#### Pathway Maps"
    # The page uses <h4> for major sections and <b> for subsections
    # Pathways are plain text nodes containing "NNNNN" followed by link text

    content = soup.find("div", {"id": "content"}) or soup.body

    for tag in content.descendants:
        # Skip non-tag nodes here; we handle NavigableString below
        if not hasattr(tag, "name"):
            continue

        # h4 → supergroup  e.g. "1. Metabolism"
        if tag.name == "h4":
            text = tag.get_text(" ", strip=True)
            # Strip leading number: "1. Metabolism" → "Metabolism"
            m = re.match(r'^\d+\.\s+(.*)', text)
            if m:
                supergroup = m.group(1).strip()
                group = ""   # reset group when supergroup changes

        # strong/b → group  e.g. "1.1 Carbohydrate metabolism"
        elif tag.name in ("b", "strong"):
            text = tag.get_text(" ", strip=True)
            m = re.match(r'^\d+\.\d+\s+(.*)', text)
            if m:
                group = m.group(1).strip()

        # a → pathway link  e.g. href contains /pathway/map00010
        elif tag.name == "a" and tag.get("href", ""):
            href  = tag["href"]
            annot = tag.get_text(" ", strip=True)

            # Match map, ko, hsa, dme, ath, sce, xla, cel, ccr, ath prefixes
            m = re.search(r'/pathway/(?:map|ko|hsa|dme|ath|sce|xla|cel|ccr|cre|eco)?'
                          r'(\d{5})', href)
            if m and annot and supergroup and group:
                pid = m.group(1)
                rows.append((pid, group, supergroup, annot))

    return rows


def parse_hierarchy_from_markdown(html):
    """
    Fallback: parse from the raw text using regex when BS4 doesn't find
    the expected tag structure (some servers return simplified HTML).
    Works on the markdown-style text fetched by web_fetch.
    """
    rows = []
    supergroup = ""
    group = ""

    # Supergroup lines: "#### 1. Metabolism" or "**1. Metabolism**"
    # Group lines:      "**1.1 Carbohydrate metabolism**"
    # Pathway lines:    "00010 M N\n[Glycolysis...](...)"

    for line in html.split("\n"):
        line = line.strip()

        # Supergroup: #### N. Name  or  **N. Name**
        m = re.match(r'^#{2,4}\s+\d+\.\s+(.+)', line)
        if m:
            supergroup = m.group(1).strip().rstrip("*").strip()
            group = ""
            continue

        # Group: **N.N Name**
        m = re.match(r'^\*{1,2}(\d+\.\d+)\s+(.+?)\*{0,2}$', line)
        if m:
            group = m.group(2).strip()
            continue

        # Pathway link: [Annotation](https://.../mapNNNNN)
        m = re.search(r'\[([^\]]+)\]\(https?://[^\)]*?(\d{5})[^\)]*\)', line)
        if m:
            annot = m.group(1).strip()
            pid   = m.group(2)
            if supergroup and group and not annot.startswith("map"):
                rows.append((pid, group, supergroup, annot))

    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", required=True,
                        help="Output TSV path")
    parser.add_argument("--url", default=URL,
                        help=f"KEGG pathway list URL (default: {URL})")
    args = parser.parse_args()

    print(f"Fetching {args.url}...", file=sys.stderr)
    html = fetch_page(args.url)

    print("Parsing pathway hierarchy...", file=sys.stderr)
    rows = parse_hierarchy(html)

    # Fallback to markdown parser if BS4 found nothing
    if not rows:
        print("  BS4 parse returned no rows — trying text parser...",
              file=sys.stderr)
        rows = parse_hierarchy_from_markdown(html)

    # Deduplicate while preserving order
    seen = set()
    deduped = []
    for row in rows:
        if row not in seen:
            seen.add(row)
            deduped.append(row)

    print(f"  {len(deduped)} pathways parsed", file=sys.stderr)

    with open(args.output, "w") as fh:
        fh.write("Pathway\tGroup\tSupergroup\tAnnot\n")
        for pid, group, supergroup, annot in deduped:
            fh.write(f"{pid}\t{group}\t{supergroup}\t{annot}\n")

    print(f"Done → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
