from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET


def executable(value):
    path = shutil.which(str(value))
    if path is None:
        raise FileNotFoundError(f"Executable not found: {value}. Supply --solver or --mpi-exec.")
    return str(Path(path).resolve())


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
    recorded_command = [Path(command[0]).name, *command[1:3], Path(command[3]).name, command[4]]
    print("Launching:", " ".join(recorded_command), flush=True)
    with (case_dir / "run.log").open("a", buffering=1) as log:
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
    # Retain run.log; discard the solver's separate iteration history only after success.
    for name in ("histor.dat", "history.dat"):
        (output / name).unlink(missing_ok=True)
    return True
