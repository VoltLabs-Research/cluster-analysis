#!/usr/bin/env python3
"""Generic Volt plugin wrapper: bridges the persistent stub protocol to a
native binary that consumes a LAMMPS dump file. STDOUT is reserved for the
stub's binary IPC protocol — never write there from Python or subprocesses."""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

BINARY_NAME = "cluster-analysis"
PLUGIN_REPO_DIRNAME = "ClusterAnalysis"
ENV_BINARY_OVERRIDE = "VOLT_CLUSTER_ANALYSIS_BINARY"
REQUIRED_OUTPUTS = ["_cluster_analysis.parquet", "_atoms.parquet"]
LOG_TAG = "cluster-analysis-plugin"

SCRIPT_DIR = Path(__file__).resolve().parent
PLUGIN_ROOT = Path(os.environ.get("PLUGIN_PROJECT_DIR", SCRIPT_DIR.parent)).resolve()
PLUGINS_ROOT = PLUGIN_ROOT.parent if PLUGIN_ROOT.parent.name == "plugins" else None
EMBEDDED_LOADER = (PLUGIN_ROOT / "lib/ld-linux-x86-64.so.2").resolve()
EMBEDDED_LIBRARY_DIR = (PLUGIN_ROOT / "lib").resolve()


class WrapperError(RuntimeError):
    pass


def _resolve_binary() -> Path:
    env_value = os.environ.get(ENV_BINARY_OVERRIDE, "").strip()
    candidates: list[Path] = []
    if env_value:
        candidates.append(Path(env_value))
    candidates.append(PLUGIN_ROOT / "bin" / BINARY_NAME)
    if PLUGINS_ROOT is not None:
        candidates.extend([
            PLUGINS_ROOT / PLUGIN_REPO_DIRNAME / "build/build/Release" / BINARY_NAME,
            PLUGINS_ROOT / PLUGIN_REPO_DIRNAME / "build-local/build/Release" / BINARY_NAME,
            PLUGINS_ROOT / PLUGIN_REPO_DIRNAME / "build-manual/build/Release" / BINARY_NAME,
        ])
    which = shutil.which(BINARY_NAME)
    if which:
        candidates.append(Path(which))
    for candidate in candidates:
        if candidate and candidate.exists():
            return candidate.resolve()
    probed = "\n".join(f"  - {c}" for c in candidates)
    raise WrapperError(f"No pude resolver binario {BINARY_NAME}. Paths probados:\n{probed}")


def _resolve_embedded_runtime_command(command: list[str]) -> list[str]:
    if not command:
        return command
    if not EMBEDDED_LOADER.exists() or not EMBEDDED_LIBRARY_DIR.exists():
        return command
    binary_path = Path(command[0]).resolve()
    try:
        binary_path.relative_to(PLUGIN_ROOT)
    except ValueError:
        return command
    if binary_path == EMBEDDED_LOADER:
        return command
    return [
        str(EMBEDDED_LOADER),
        "--library-path", str(EMBEDDED_LIBRARY_DIR),
        str(binary_path),
        *command[1:],
    ]


def _run(command: list[str]) -> None:
    command = _resolve_embedded_runtime_command(command)
    sys.stderr.write(f"[{LOG_TAG}] {' '.join(shlex.quote(part) for part in command)}\n")
    sys.stderr.flush()
    completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if completed.returncode != 0:
        raise WrapperError(f"El comando fallo con exit code {completed.returncode}: {command[0]}")


def _require_outputs(output_base: str) -> None:
    for suffix in REQUIRED_OUTPUTS:
        expected = Path(f"{output_base}{suffix}")
        if not expected.exists():
            raise WrapperError(f"Falta el archivo requerido: {expected}")


_VOLT_RUNTIME_FLAGS_WITH_VALUE = {
    "--selectedTimesteps",
    "--selected-timesteps",
}


def _filter_runtime_flags(args: list[str]) -> list[str]:
    filtered: list[str] = []
    i = 0
    while i < len(args):
        token = str(args[i])
        if token in _VOLT_RUNTIME_FLAGS_WITH_VALUE:
            i += 2
            continue
        filtered.append(token)
        i += 1
    return filtered


def _ensure_executable(path: Path) -> None:
    try:
        mode = path.stat().st_mode
        if not (mode & 0o111):
            path.chmod(mode | 0o755)
    except OSError:
        pass


def _run_binary_with_args(args: list[str]) -> dict:
    args = _filter_runtime_flags([str(a) for a in args])
    if len(args) < 2:
        raise WrapperError("Se esperaban al menos 2 argumentos: <input_dump> <output_base>")
    output_base = args[1]
    Path(output_base).parent.mkdir(parents=True, exist_ok=True)
    binary = _resolve_binary()
    _ensure_executable(binary)
    if EMBEDDED_LOADER.exists():
        _ensure_executable(EMBEDDED_LOADER)
    command = [str(binary), *args]
    _run(command)
    _require_outputs(output_base)
    return {
        "ok": True,
        "outputBase": output_base,
        "binary": str(binary),
        "outputs": [f"{output_base}{suffix}" for suffix in REQUIRED_OUTPUTS],
    }


def process(frame, config):
    """Volt persistent plugin entrypoint."""
    del frame
    if not isinstance(config, dict):
        raise WrapperError("config debe ser un dict")
    args = config.get("args")
    if not isinstance(args, list):
        raise WrapperError("config['args'] debe ser una lista de strings")
    return _run_binary_with_args(args)


def _main_cli() -> int:
    argv = sys.argv[1:]
    _run_binary_with_args(argv)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(_main_cli())
    except WrapperError as error:
        sys.stderr.write(f"[{LOG_TAG}] error: {error}\n")
        raise SystemExit(1)
