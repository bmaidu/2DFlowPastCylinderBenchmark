from pathlib import Path
import xml.etree.ElementTree as ET

TEMPLATE = Path(__file__).with_name("solver_template.xml")


def create_solver(case_dir, config):
    root = ET.parse(TEMPLATE).getroot()
    s = config["solver"]
    values = {
        "GeneralSimulationParameters/Continue_previous_simulation": "false",
        "GeneralSimulationParameters/Save_results_in_folder": "solver_output",
        "GeneralSimulationParameters/Time_step_size": s["dt"],
        "GeneralSimulationParameters/Number_of_time_steps": s["n_steps"],
        "GeneralSimulationParameters/Increment_in_saving_VTK_files": s["vtk_increment"],
        "GeneralSimulationParameters/Increment_in_saving_restart_files": s["restart_increment"],
        "Add_equation/Density": s["density"],
        "Add_equation/Viscosity/Value": s["viscosity"],
        "Add_equation/Add_BC[@name='inlet']/Value": -s["speed"],
    }
    for xpath, value in values.items():
        node = root.find(xpath)
        if node is None:
            raise ValueError(f"Solver template is missing {xpath}.")
        node.text = f" {value} "
    path = Path(case_dir) / "solver.xml"
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)
    return path
