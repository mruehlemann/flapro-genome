from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Iterable


LOG = logging.getLogger("flapro_genome")


class FlaProGenomeError(RuntimeError):
    """Expected user-facing error."""


def configure_logging(verbose: bool = False, quiet: bool = False) -> None:
    level = logging.WARNING if quiet else logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s", force=True)


def require_executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise FlaProGenomeError(
            f"Required executable '{name}' was not found on PATH. "
            "Install the flapro-genome Conda environment first."
        )
    return executable


def run_command(command: Iterable[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    args = [str(item) for item in command]
    LOG.debug("Running: %s", " ".join(args))
    try:
        return subprocess.run(
            args,
            cwd=cwd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except FileNotFoundError as exc:
        raise FlaProGenomeError(f"Executable not found: {args[0]}") from exc
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or "").strip()
        raise FlaProGenomeError(
            f"Command failed with exit code {exc.returncode}: {' '.join(args)}"
            + (f"\n{details}" if details else "")
        ) from exc


def command_version(command: list[str]) -> str:
    try:
        result = run_command(command)
    except FlaProGenomeError:
        return "unknown"
    text = (result.stdout or result.stderr).strip()
    return text.splitlines()[0] if text else "unknown"


def safe_stem(path: Path) -> str:
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in (".fasta", ".faa", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem

