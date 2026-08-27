"""Download + raw-file plumbing shared by all loaders."""
from __future__ import annotations

import hashlib
from pathlib import Path

import requests


class RawMissingError(FileNotFoundError):
    """Raised when a required raw file is absent, with a fix hint."""

    def __init__(self, path: Path, source: str, how: str):
        super().__init__(
            f"Missing raw file for {source}: {path}\n"
            f"  -> {how}\n"
            f"  (or run the synthetic demo: `warhead demo`)"
        )


def sha256_of(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download_file(
    url: str,
    dest: str | Path,
    *,
    sha256: str | None = None,
    overwrite: bool = False,
    timeout: int = 60,
) -> Path:
    """Stream ``url`` to ``dest``, verify checksum if given, skip if present.

    Raw sources are checksummed on download (WARHEAD.md sec 4): a mismatch raises
    rather than silently using a corrupt/rotated file.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not overwrite:
        if sha256 and sha256_of(dest) != sha256:
            raise ValueError(f"{dest} exists but sha256 mismatch; delete to re-download")
        return dest

    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        tmp.replace(dest)

    if sha256:
        got = sha256_of(dest)
        if got != sha256:
            raise ValueError(f"sha256 mismatch for {dest}: expected {sha256}, got {got}")
    return dest
