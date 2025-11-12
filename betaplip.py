"""
PyMOL plugin with functionalities from PLIP

    Imstepf, Nicolas. Re-Programming of Fe(II)/a-Ketoglutarate Dependent Hydroxylases to Halogenases
    Master Thesis (Zurich University of Applied Sciences, 2023).

    Adasme, M. F. et al. PLIP 2021: expanding the scope of the protein-ligand 
    interaction profiler to DNA and RNA. Nucleic Acids Res. 49, W530-W534 (W1 2021).

betaplip.py (equivalent with main.py)
"""


__author__      = "David Patsch & Nicolas Imstepf"
__copyright__   = "Copyright 2023, Zurich University of Applied Sciences, Competence Center for Biocatalysis"
__credits__     = ["Rebecca Buller", "Sean Hueppi", "Fabian Meyer", "Michael Niklaus"]
__license__     = "GPL3, images: CC BY-NC-SA 4.0"
__version__     = "1.0.0"
__maintainer__  = "Nicolas Imstepf"
__email__       = "nicolasimstepf@gmail.com"


import os
import sys

# add betplip path to directory 
sys.path.insert(1, os.path.dirname(__file__))
from component import InteractionGUI


def __init_plugin__(app=None) -> None:
    """
    helper function to add betaplip to PyMOL menu
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('BetaPlip', run_plugin_gui)

# global reference to avoid garbage collection of our dialog
dialog = None

def run_plugin_gui() -> None:
    """
    helper function to run the GUI
    """
    global dialog
    if dialog is None:
        dialog = InteractionGUI()
    dialog.show()

# show application window by running betaplip.py (developing purpose)
if __name__ == "__main__":
    from PyQt5.QtWidgets import QApplication
    print("MAIN")
    app = QApplication(sys.argv)
    wnd = InteractionGUI()
    wnd.show()
    sys.exit(app.exec_())
