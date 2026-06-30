# ingest.py — PDF -> Markdown ingestion for the Chemdem drop-zone panel.
# Wraps the same pymupdf4llm conversion as data/ingest_pdf.py, but operates on
# uploaded bytes and saves into data/papers/.

import os, re, datetime

_BACKEND = os.path.dirname(os.path.abspath(__file__))
# Writable papers dir. Override with CHEMDEM_PAPERS_DIR on hosted deploys
# (e.g. Hugging Face Spaces) where the repo-relative ../data path may not exist.
PAPERS = os.environ.get("CHEMDEM_PAPERS_DIR") or os.path.abspath(
    os.path.join(_BACKEND, "..", "data", "papers"))


def _safe_stem(filename: str) -> str:
    stem = os.path.splitext(os.path.basename(filename or "paper"))[0]
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_")
    return stem or "paper"


def convert_pdf_bytes(filename: str, data: bytes) -> dict:
    """
    Save the uploaded PDF to data/papers/, convert it to Markdown, save the .md
    alongside it, and return metadata + a short preview.
    """
    import pymupdf4llm

    os.makedirs(PAPERS, exist_ok=True)
    stem = _safe_stem(filename)
    pdf_path = os.path.join(PAPERS, stem + ".pdf")
    md_path = os.path.join(PAPERS, stem + ".md")

    with open(pdf_path, "wb") as f:
        f.write(data)

    md = pymupdf4llm.to_markdown(pdf_path)
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md)

    words = len(md.split())
    table_rows = md.count("|---")
    preview = md.strip()[:1200]

    return {
        "ok": True,
        "pdf_filename": os.path.basename(pdf_path),
        "md_filename": os.path.basename(md_path),
        "words": words,
        "table_rows": table_rows,
        "chars": len(md),
        "preview": preview,
        "saved_at": datetime.datetime.now().isoformat(timespec="seconds"),
    }


def _title(text: str) -> str:
    for line in text.splitlines():
        s = line.strip().lstrip("#").strip().strip("*").strip()
        if len(s) > 3:
            return s[:140]
    return "Untitled paper"


def _abstract(text: str) -> str:
    # first chunk of prose after the title, markdown markers lightly stripped
    lines = [l.strip() for l in text.splitlines()]
    body = []
    seen_title = False
    for l in lines:
        if not l:
            continue
        if not seen_title:
            seen_title = True
            continue
        if l.startswith("|") or l.startswith("==>") or l.startswith("!["):
            continue
        body.append(l.lstrip("#").strip())
        if sum(len(x) for x in body) > 360:
            break
    return " ".join(body)[:380]


def list_papers() -> list:
    """List converted markdown files already in the papers dir, with previews."""
    if not os.path.isdir(PAPERS):
        return []
    out = []
    for name in sorted(os.listdir(PAPERS)):
        if not name.lower().endswith(".md") or name.lower() == "readme.md":
            continue
        path = os.path.join(PAPERS, name)
        try:
            text = open(path, encoding="utf-8").read()
        except Exception:
            text = ""
        out.append({
            "md_filename": name,
            "title": _title(text),
            "abstract": _abstract(text),
            "words": len(text.split()),
            "table_rows": text.count("|---"),
            "has_pdf": os.path.exists(os.path.join(PAPERS, name[:-3] + ".pdf")),
            "modified": datetime.datetime.fromtimestamp(
                os.path.getmtime(path)).isoformat(timespec="seconds"),
        })
    return out


def file_path(md_filename: str) -> str | None:
    """Resolve a download path for a paper: the original PDF if kept, else the .md."""
    name = os.path.basename(md_filename or "")
    stem = name[:-3] if name.lower().endswith(".md") else name
    pdf = os.path.join(PAPERS, stem + ".pdf")
    md = os.path.join(PAPERS, stem + ".md")
    if os.path.isfile(pdf):
        return pdf
    if os.path.isfile(md):
        return md
    return None


def get_markdown(md_filename: str) -> str | None:
    """Return the full markdown text of a converted paper, or None if missing."""
    name = os.path.basename(md_filename or "")
    if not name.lower().endswith(".md"):
        name += ".md"
    path = os.path.join(PAPERS, name)
    if not os.path.isfile(path):
        return None
    try:
        return open(path, encoding="utf-8").read()
    except Exception:
        return None
