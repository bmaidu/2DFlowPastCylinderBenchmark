from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime, timezone


def executable(value):
    path = shutil.which(str(value))
    if path is None:
        raise FileNotFoundError(f"Executable not found: {value}. Supply --solver or --mpi-exec.")
    return str(Path(path).resolve())


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def final_step_reached(case_dir):
    """Check only whether a nonempty VTU snapshot reaches the XML's final step."""
    case_dir = Path(case_dir)
    solver_xml = case_dir / "solver.xml"
    if not solver_xml.is_file():
        return False
    general = ET.parse(solver_xml).getroot().find("GeneralSimulationParameters")
    final_step = int(general.findtext("Number_of_time_steps"))
    output = case_dir / general.findtext("Save_results_in_folder").strip()
    prefix = general.findtext("Name_prefix_of_saved_VTK_files").strip()
    for path in output.glob(f"{prefix}_*.vtu"):
        step = path.stem[len(prefix) + 1:]
        if (step.isdigit() and int(step) >= final_step
                and path.is_file() and path.stat().st_size > 0):
            return True
    return False


def run_solver(case_dir, solver="svmultiphysics", mpi_exec="mpirun", np=8):
    case_dir = Path(case_dir).resolve()
    if np < 1:
        raise ValueError("--np must be positive.")
    if final_step_reached(case_dir):
        print("Final VTU timestep has been reached; skipping the solver.", flush=True)
        return False
    stop = case_dir / "STOP_SIM"
    if stop.exists():
        raise FileExistsError("STOP_SIM exists; remove it explicitly before starting a fresh run.")
    command = [executable(mpi_exec), "-n", str(np), executable(solver), "solver.xml"]
    output = case_dir / ET.parse(case_dir / "solver.xml").getroot().findtext(
        "GeneralSimulationParameters/Save_results_in_folder").strip()
    output.mkdir(parents=True, exist_ok=True)
    provenance_path = case_dir / "provenance.json"
    provenance = json.loads(provenance_path.read_text()) if provenance_path.exists() else {}
    recorded_command = [Path(command[0]).name, *command[1:3], Path(command[3]).name, command[4]]
    provenance["execution"] = {
        "command": recorded_command, "cwd": ".",
        "solver_sha256": sha256(command[3]),
        "started_utc": datetime.now(timezone.utc).isoformat(), "status": "running",
    }
    provenance_path.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    print("Launching:", " ".join(recorded_command), flush=True)
    with (case_dir / "run.log").open("a", buffering=1) as log:
        try:
            with subprocess.Popen(command, cwd=case_dir, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True, errors="replace",
                                  bufsize=1) as process:
                try:
                    for line in process.stdout:
                        log.write(line)
                        print(line, end="", flush=True)
                    returncode = process.wait()
                except BaseException:
                    process.kill()
                    process.wait()
                    raise
                if returncode:
                    raise subprocess.CalledProcessError(returncode, command)
        except BaseException as error:
            reason = str(error).replace(str(case_dir), ".")
            for index in (0, 3):
                reason = reason.replace(command[index], Path(command[index]).name)
            provenance["execution"].update(status="failed", error=reason,
                finished_utc=datetime.now(timezone.utc).isoformat())
            provenance_path.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
            raise
    provenance["execution"].update(status="completed", returncode=0,
        finished_utc=datetime.now(timezone.utc).isoformat())
    provenance_path.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return True
