"""
Graphical user interface

betaplipgui.py
"""

from . import api, styles
from .components import ComboBox, MyBar


from pathlib import Path
from functools import partial
import re
import os
import shutil

from PyQt5.QtWidgets import (QDialog, QVBoxLayout,QHBoxLayout, QPushButton, 
                            QGroupBox, QGridLayout, QLabel, QCheckBox, QLineEdit, 
                            QFileDialog, QSpinBox, QWidget, QRadioButton)
from PyQt5.QtCore import Qt, QDir, QTimer

from pymol import cmd


class InteractionGUI(QDialog):
    """Plugin to run Plip & show interactions"""

    def __init__(self, parent=None) -> None:
        super(InteractionGUI, self).__init__(parent)

        # Multithreading: speed up pymol in general (performing geometry builds in parallel)
        cmd.set("async_builds", 1)      

        # initialize list of pymol objects (molecules only) & API
        self.selections = self.getMolecules()
        # self.API = api.plipAPI(Path(r"C:\Users\imst\Temp")) #TODO: bugfixing  nothing is working anymore 250919
        self.API = api.plipAPI(Path(""))

        # add icon image path
        # from stackoverflow.com/questions/70624578
        root = os.path.dirname(os.path.abspath(__file__))
        QDir.addSearchPath('images', os.path.join(root, 'images'))

        # Layout
        self.setStyleSheet(styles.David)
        self.setMinimumSize(500, 800)
        # self.setFixedSize(600, 850)
        self.setWindowFlags(Qt.FramelessWindowHint)

        # Main Layout
        MainLayout = QVBoxLayout()
        self.MenuBar = MyBar(self)
        self.MenuBar.actionshowDist.triggered.connect(self.showDist)
        self.MenuBar.actionclean_up.triggered.connect(self.clean_up)
        MainLayout.addWidget(self.MenuBar)
        MainLayout.setContentsMargins(10,10,10,10)
        
        # Input
        InputGroupBox = QGroupBox("Input", objectName="Input")
        InputLayout = QGridLayout()

        self.alphaDockBtn = QRadioButton("AlphaDock")
        self.alphaDockBtn.setChecked(True)
        
        self.pdbBtn = QRadioButton("PDB files")
        self.pdbBtn.toggled.connect(self.setPDBonly)

        # Input from Pymol Object
        self.pymolBtn = QRadioButton("Pymol Object")
        self.inputComboBox = ComboBox(objectName="InputObjects")
        self.inputComboBox.addItems(self.selections)
        self.inputComboBox.entered.connect(self.updateInput)
        self.inputComboBox.setEnabled(False)
        self.pymolBtn.toggled.connect(lambda: self.inputComboBox.setEnabled(self.pymolBtn.isChecked()))

        # Load Directory Button -> open QFileDialog
        loadDirectoryBtn = QPushButton("Load")
        loadDirectoryBtn.setStyleSheet(styles.DirectoryButton)
        loadDirectoryBtn.setToolTip(
                        "Load the directory containing the pdb files.\n"
                        "For example result folder '1' from AlphaFold.\n"
                        "Do not load the 'plip' folder directly." )
        loadDirectoryBtn.clicked.connect(self.setDirectory)

        # Display path
        self.DirectoryLabel = QLabel()

        # Plip Button -> run Plip
        analyzeDirectoryBtn = QPushButton("Plip!")
        analyzeDirectoryBtn.setToolTip(
                            "Analyse the result with Plip. Load the\n" 
                            "resulting interactions into PyMOL.")
        analyzeDirectoryBtn.setStyleSheet(styles.HoverGreenButton)
        analyzeDirectoryBtn.clicked.connect(self.runPlip)

        InputLayout.addWidget(loadDirectoryBtn, 0, 0)
        InputLayout.addWidget(self.DirectoryLabel, 0, 1, 1, 2)
        InputLayout.addWidget(self.alphaDockBtn, 0, 3)
        InputLayout.addWidget(self.pdbBtn, 1, 3)
        InputLayout.addWidget(self.pymolBtn, 3, 0)
        InputLayout.addWidget(analyzeDirectoryBtn, 3, 3)
        InputLayout.addWidget(self.inputComboBox, 3, 1)
        InputGroupBox.setLayout(InputLayout)
        MainLayout.addWidget(InputGroupBox)

        # Distance
        DistanceGroupBox = QGroupBox("Distance Filter")
        DistanceLayout = QGridLayout()
        
        # Selections (Pymol Objects: Molecules)
        for i in range(2):
            selectionLabel = QLabel(f"Selection {i+1}") 
            selComboBox = ComboBox(objectName=f"Sel_{i+1}")
            selComboBox.addItems(self.selections)
            selComboBox.entered.connect(self.updateSelections)
            DistanceLayout.addWidget(selectionLabel, i+1, 0)
            DistanceLayout.addWidget(selComboBox, i+1, 1)

        # Distance Filter
        DistanceLabel = QLabel("Distance (\u00C5)")
        self.DistanceLEdit = QLineEdit()
        self.DistanceLEdit.editingFinished.connect(self.setDistFilter)
        DistanceLayout.addWidget(DistanceLabel, i+2, 0)
        DistanceLayout.addWidget(self.DistanceLEdit, i+2, 1)
        DistanceGroupBox.setLayout(DistanceLayout)
        MainLayout.addWidget(DistanceGroupBox)

        # Interactions
        InteractionGroupBox = QGroupBox("Interactions")
        InteractionLayout = QGridLayout()

        # Description
        btnSelectAll = QCheckBox("Select/Deselect all", objectName="selectAll")
        btnSelectAll.setChecked(True)
        btnSelectAll.clicked.connect(self.selectAll)
        InteractionLayout.addWidget(btnSelectAll)
        InteractionLayout.addWidget(QLabel('Residue (e.g., 1, 2, 3)'), *(0, 1))

        # Choice / Interaction Types
        # Debouncing for the checkboxes
        # adapted from gist.github.com/chipolux/a600d2a31b6811d553651822f89c9e39
        self.debounceTypes = QTimer()
        self.debounceTypes.setInterval(500)
        self.debounceTypes.setSingleShot(True)
        self.debounceTypes.timeout.connect(self.setInteractionFilter)

        self.debounceMinAnd = QTimer()
        self.debounceMinAnd.setInterval(500)
        self.debounceMinAnd.setSingleShot(True)
        self.debounceMinAnd.timeout.connect(self.setMinAndInt)

        for i, (k, v) in enumerate(styles.PYMOL_STYLES.items()):
            btn = QCheckBox(k, objectName=k)
            btn.setToolTip(v["description"])
            btn.setChecked(True)
            self.API.interactionFilter[k] = True
            btn.toggled.connect(self.debounceTypes.start)
            InteractionLayout.addWidget(btn, *(i+1, 0))

            ledit = QLineEdit(objectName=f"{k}_res")
            ledit.editingFinished.connect(partial(self.setResidueFilter, k))
            InteractionLayout.addWidget(ledit, *(i+1, 1))

        # Logical Filter Buttons
        layoutFilterBtn = QHBoxLayout()

        # OR
        self.buttonOR = QPushButton("OR")
        self.buttonOR.clicked.connect(partial(self.setLogicalFilter, "OR"))
        self.buttonOR.setStyleSheet(styles.LogicalButton)
        layoutFilterBtn.addWidget(self.buttonOR)

        # AND
        self.buttonAND = QPushButton("AND")
        self.buttonAND.clicked.connect(partial(self.setLogicalFilter, "AND"))
        layoutFilterBtn.addWidget(self.buttonAND)

        self.logicalPushButtons = {
            "OR":  self.buttonOR,
            "AND": self.buttonAND
        }

        InteractionLayout.addLayout(layoutFilterBtn, i+2, 0, 1, 4)
        InteractionGroupBox.setLayout(InteractionLayout)

        # Show interactions for selection of residues
        intraResLabel = QLabel("Highlight Residues")
        self.intraResidueFilter = None
        self.intraResidueLEdit = QLineEdit()
        self.intraResidueLEdit.editingFinished.connect(self.setIntraRes)
        InteractionLayout.addWidget(intraResLabel, i+3, 0)
        InteractionLayout.addWidget(self.intraResidueLEdit, i+3, 1)

        # hide empty checkbox
        self.hideEmptyCheckBox = QCheckBox("hide states without interactions")
        self.hideEmptyCheckBox.setChecked(self.API.removeEmpty)
        self.hideEmptyCheckBox.toggled.connect(self.setHideEmpty)
        
        layoutHideAndBtns = QHBoxLayout()
        layoutHideAndBtns.addWidget(self.hideEmptyCheckBox)

        # minimum interactions for AND operation
        self.minAndIntCheckBox = QCheckBox("min No.  of types  ")
        self.minAndIntCheckBox.setToolTip(
                                        "Set a minimum number of different interaction types.\n"
                                        "For example 3 out of 4 selected types must occur.")
        self.minAndIntCheckBox.toggled.connect(self.activateSpinBox)
        self.minAndIntCheckBox.setEnabled(False)
        self.minAndInteractions = QSpinBox()
        self.minAndInteractions.setSpecialValueText(" ")
        self.minAndInteractions.valueChanged.connect(self.debounceMinAnd.start)
        self.minAndInteractions.setEnabled(False)

        layoutAnd = QHBoxLayout()
        layoutAnd.addWidget(self.minAndIntCheckBox)
        layoutAnd.addWidget(self.minAndInteractions)
        layoutHideAndBtns.addLayout(layoutAnd)

        InteractionLayout.addLayout(layoutHideAndBtns, i+4, 0, 1, 4)

        InteractionGroupBox.setLayout(InteractionLayout)
        MainLayout.addWidget(InteractionGroupBox)

        # Outputs
        OutputGroupBox = QGroupBox("Outputs")
        OutputLayout = QGridLayout()

        # load button
        loadBtn = QPushButton("Load Filter") 
        loadBtn.clicked.connect(self.loadFilter)
        OutputLayout.addWidget(loadBtn, 0, 0)

        # filename
        self.filtflnmLEdit = QLineEdit()
        self.filtflnmLEdit.setPlaceholderText("filter_filename")
        OutputLayout.addWidget(self.filtflnmLEdit, 0, 1, 1, 2)

        # save button
        saveBtn = QPushButton("Save Filter")
        saveBtn.clicked.connect(self.saveFilter)
        saveBtn.setStyleSheet(styles.HoverPurpleButton)
        OutputLayout.addWidget(saveBtn, 0, 3)

        OutputGroupBox.setLayout(OutputLayout)
        MainLayout.addWidget(OutputGroupBox)
        self.setLayout(MainLayout)

        # disable GroupBoxes as no directory is set
        for box in self.findChildren(QGroupBox):
            if box.objectName() != "Input":
                box.setEnabled(False)
        
        self.enableEverything(False)

        self.setWindowTitle("Pliper")
        self.show()
    
    def enableEverything(self, enable: bool=True) -> None:
        """
        enables all group boxes, when result exists
        exception: input groupbox, that is always 
        """
        for groupbox in self.findChildren(QGroupBox):
            if groupbox.objectName() != "Input":
                groupbox.setEnabled(enable)

    def change_host(self, new_host: str, configs: dict) -> None:
        """
        changes host
        called by QAMenuHosts.change_host_and_update_ticks() from components.py

        :param str new_host:                        hostname
        :param dict[str[Union[str, int]]] configs:  host parameters (alias, num_cpu, port)
        """
        self.API.hostname = new_host

        if "num_cpu" in configs.keys():
            self.API.cpu_cnt = configs["num_cpu"]
        else:
            self.API.cpu_cnt = 4

        if "connect_on_port" in configs.keys():
            self.API.connect_on_port = configs["connect_on_port"]
        else:
            self.API.connect_on_port = 77

        if "username" in configs.keys():
            self.API.username= configs["username"]

        if "password" in configs.keys():
            self.API.password = configs["password"]

        print("Host changed to:", self.API.hostname)
        print("Port changed to:", self.API.connect_on_port)


    @staticmethod
    def str2float(value: str) -> float:
        """
        transform an input string into a floating point

        exception handling for distance input 
        """
        # exception handling if value cannot be transformed
        if value:
            try:
                value = float(value)
                # additional information if value is negative
                if value <= 0:
                    print(f"Warning: Distance value is {value}")
                return value
            except:
                print("Distance must be a number")

    def setIntraRes(self) -> None:
        """
        set intraRes based on input
        """
        inputtxt = self.intraResidueLEdit.text()
        if inputtxt == "":
            self.intraResidueFilter = None
        else:
            residues = inputtxt.split(",")
            self.intraResidueFilter = [str(self.str2int(x)) for x in residues]
        self.showInteractions()

    def setDistFilter(self) -> None:
        """
        set distFilter based on Selection1/2 and distance
        """
        # reset distFilter and initialize placeholder
        self.API.distFilter = {}
        selection = {}

        # get atom coordinates from selections
        for i in range(2):
            # get selected pymol objects
            cbox = self.findChild(ComboBox, f"Sel_{i+1}")
            pymolObject = cbox.currentText()

            # iterate over different states and get selected atom coordinates
            if pymolObject != None and pymolObject:
                nstates = cmd.count_states(pymolObject)
                for state in range(1, nstates+1):
                    atoms = cmd.get_model(pymolObject, state)
                    filename = cmd.get_title(pymolObject, state)
                    if atoms.atom:
                        selection.setdefault(f"Sel_{i+1}", {})
                        selection[f"Sel_{i+1}"][filename] = [(at.coord[0], at.coord[1], at.coord[2]) for at in atoms.atom]
        
        # get entered distance
        distance = self.str2float(self.DistanceLEdit.text())
        if distance:
            selection["distance"] = distance

        # if selection is valid, overwrite distFilter and show interactions
        if len(selection) == 3:
            self.API.distFilter = selection
            self.showInteractions()
            
            cmd.distance("Distance_Filter", self.findChild(ComboBox, f"Sel_1").currentText(), self.findChild(ComboBox, f"Sel_2").currentText(), cutoff = selection["distance"])
            
    def updateSelections(self) -> None:
        """
        update Selection comboboxes with pymol objects (molecules only)
        """
        for i in range(2):
            self.selections = self.getMolecules()
            cbox = self.findChild(ComboBox, f"Sel_{i+1}")
            item = cbox.currentText()
            cbox.clear()       
            cbox.addItems(self.selections)
            cbox.setCurrentText(item)

    def updateInput(self) -> None:
        """
        update Input combobox with pymol objects (molecules only)
        """
        self.selections = self.getMolecules()
        item = self.inputComboBox.currentText()
        self.inputComboBox.clear()       
        self.inputComboBox.addItems(self.selections)
        self.inputComboBox.setCurrentText(item)

    @staticmethod
    def getMolecules() -> list:
        """
        get all pymol objects of class molecule
        :return:    list of strings
        """
        return ["None"] + [pymolObject for pymolObject in cmd.get_names("all") if cmd.get_type(pymolObject) == "object:molecule"]

    def setDirectory(self) -> None:
        """
        open FileDialog to set directory
        call self.API.updateDirectory()
        """
        # select directory via QFileDialog
        options = QFileDialog.Options()
        directory = QFileDialog.getExistingDirectory(self, "select the directory",
                                                     options=options) + "/"

        if directory == "/":
            return

        # show path next to the button
        # exception handling / shorten long paths
        idx = 0
        if len(directory) > 30:
            # adapted from stackoverflow.com/questions/4664850
            slashpos = [m.start() for m in re.finditer('/', directory)]
            for i in reversed(slashpos):
                if i < len(directory) - 30:
                    idx = i
                    break

        self.DirectoryLabel.setText("..." + directory[idx:])
        self.API.updateDirectory(Path(directory))

    def runPlip(self) -> None:
        """
        call self.API.run_plip() and show interactions
        """  
        # exception handling if no directory is set
        if self.API.directory == Path("") and (self.pdbBtn.isChecked() or self.alphaDockBtn.isChecked()):
            print("Please select a directory")
            return
        elif self.pymolBtn.isChecked() and self.inputComboBox.currentText() == "None":
            print("Please select a Pymol Object")
            return

        # exception handling if Structure is selected from PyMOL
        if self.pymolBtn.isChecked() and self.inputComboBox.currentText() != "None":
            if self.API.directory == Path("") or self.API.directory == None:
                tempfolder = Path.home() / "temp_Plip"
                tempfolder.mkdir(exist_ok=True)
                print(tempfolder, "created")
                self.API.updateDirectory(tempfolder)

            selection = self.inputComboBox.currentText()
            cmd.save(self.API.directory / "temp_selection.pdb", selection)
            self.pymolOutProcessing(self.API.directory / "temp_selection.pdb")
            self.API.pdbonly = True
        
        # run Plip
        self.API.run_plip()
        
        # remove temporary folder
        if self.API.directory.name == "temp_Plip":
            shutil.rmtree(self.API.directory)
            print(tempfolder, "removed")
            self.API.directory = None
            self.API.plipfolder = None
            self.API.resultfile = None

        if self.API.results:
            self.enableEverything()
            self.showInteractions()
            cmd.zoom("results")

    @staticmethod
    def list_rindex(li: list, x: str) -> int:
        """
        similar to str.rfind(x)
        adapted from https://stackoverflow.com/questions/6890170/

        Could be solved as one-liner
            next(i for i in reversed(range(len(li))) if li[i] == 'a')
        """
        for i in reversed(range(len(li))):
            if li[i].startswith(x):
                return i
        raise ValueError("{} is not in list".format(x))

    def pymolOutProcessing(self, pdbfile: str) -> None:
        """
        changes order of objects in pdb file
        (pymol saves HETATM's chainwise, can cause problems during plip analysis)
        """
        # get file content
        with open(pdbfile, "r") as file:
            lines = file.readlines()

        # split content ("remove" hetatm's)
        # adapted from https://stackoverflow.com/questions/949098/
        miscellaneous, hetatm = [], []
        for line in lines:
            (miscellaneous, hetatm)[line.startswith("HETATM")].append(line)

        # insert hetatm's after last TER
        ter_idx = self.list_rindex(miscellaneous, "TER")

        # overwrite file
        with open(pdbfile, "w") as file:
            for idx, line in enumerate(miscellaneous):
                # insert hetatm's after last TER
                if idx == ter_idx + 1:
                    file.write("".join(hetatm))
                file.write(line)

    def setPDBonly(self) -> None:
        """
        if only pdb files should be analyzed
        """
        self.API.pdbonly = self.pdbBtn.isChecked()

    def style(self, interactions: dict) -> None:
        """
        styling according the different interaction types
        """
        cmd.show_as("sticks", "receptor_residues")   
        self.ngl(selection="receptor_residues")

        for k, v in styles.PYMOL_STYLES.items():
            if k not in interactions:
                continue
            for sub_k, sub_v in v["style"].items():
                cmd.set(sub_k, sub_v, k)

    @staticmethod
    def cheat_interactions(interactions: list, num_states: int) -> None:
        """
        exception handling if not all states contain interactions

        :param list[str]:       list of all different interaction types
        :param int num_states:  number of states
        """
        for int in interactions:
            for i in range(num_states):
                cmd.pseudoatom("_pseudos1",
                               pos=(["100000", "100000", "100000"]), state=i+1)

                cmd.pseudoatom(f"_pseudos2",
                               pos=(["100000", "100000", "100001"]), state=i+1)

                cmd.distance(int, "_pseudos1", "_pseudos2", state=i+1)
                cmd.delete("_pseudos1, _pseudos2")

    @staticmethod
    def filterCoords(inter: dict, inner_key: str) -> tuple:
        """
        return protcoo & ligcoo values from single interaction dictionary

        :param dict[str, str] inter:    single interaction from results dictionary (among others: containing coordinates)
        :param str inner_key:           interactiontype e.g., "metal_complex"
        :return:                        return tuple of dicts (prot_coords, lig_coords)
        """
        # exception handling for metal complexes
        if inner_key != "metal_complex":
            return inter["protcoo"], inter["ligcoo"]
        return inter["metalcoo"], inter["targetcoo"]

    @staticmethod
    def getCoords(coord: dict) -> list:
        """
        get coordinates out of a dictionary

        :param dict[str: str] coord:    coordinates e.g., {"x": "0", "y": 0, "z": 0}
        :return:                        list[float] e.g., [0.0, 0.0, 0.0]
        """
        return [float(coord["x"]),
                float(coord["y"]),
                float(coord["z"])]

    def showInteractions(self) -> None:
        """
        show interactions of filtered results
        """   
        # initialize view and pymol objects
        myview = cmd.get_view()
        cmd.undo_disable()
        cmd.delete("plip")
        cmd.disable("all")

        # fix state (if not, changing this lead to bugs)
        cmd.set("state", 1)

        # set color
        cmd.set("cartoon_color","white")
        cmd.set("stick_color", "atomic")

        # filter 
        filtered_dict = self.API.filterResults()

        # exception handling if no structure fullfill criteria
        if filtered_dict == {}:
            cmd.set_view(myview)
            cmd.undo_enable()
            return

        # initialize list of unique_interactions (different types for styling)
        unique_interactions = []
        # initialize number of interactions
        all_interactions = list(styles.PYMOL_STYLES.keys())
        total_num_models = len(filtered_dict.keys())

        interactions_to_keep = [x for x in self.API.interactions if self.API.interactionFilter[x]]

        # exception handling if not all states contain interactions
        self.cheat_interactions(all_interactions, total_num_models)

        # iterate over different models & their interactions
        for t, v in enumerate(filtered_dict.values()):
            # finish & discrete should increase efficiency 
            # pymol.org/dokuwiki/doku.php?id=command:load
            if self.pymolBtn.isChecked():
                cmd.create(f"plip{t}", self.inputComboBox.currentText())
            else:
                cmd.load(self.API.plipfolder / v["path"], f"plip{t}", finish=1, discrete=1)
        
            # color ligands/ other hetatms
            cmd.color("green", "hetatm and elem C")

            interactions = v["interactions"].copy()

            # exception handling for receptor_residues (have a selection in every state)
            cmd.pseudoatom("tempSeed", pos=[0, 0, 0])
            cmd.select("temp01927824", "tempSeed", merge=True)

            # iterate over different interaction types

            for k, v_inner in interactions.items():
                # k is interaction type, eg hydrophobic_interactions
                # v is either dict with interaction or list of dicts

                # exception handling (no interaction present / interacton should be filtered)
                if not v_inner or k not in interactions_to_keep:
                    continue

                unique_interactions.append(k)

                # iterate over different interaction types
                for inner_key, vals in v_inner.items():
                    # vals are [dict] or dict
                    if type(vals) == dict:
                        vals = [vals]

                    # iterate over different interactions of a single type
                    for inter in vals:
                        if type(inter) != dict:
                            continue

                        # select all residues on receptor for receptor_residues object
                        
                        if self.intraResidueFilter \
                                and (inter["resnr_lig"] not in self.intraResidueFilter \
                                and inter["resnr"] not in self.intraResidueFilter):
                                    continue                        
                        try:
                            cmd.select(
                                f"temp01927824", f"resid {inter['resnr']} and chain {inter['reschain']} and model plip{t}", merge=True)

                            cmd.select(
                                f"temp01927824", f"resid {inter['resnr_lig']} and chain {inter['reschain_lig']} and model plip{t}", merge=True)
                        
                        except:
                            pass

                        # get protein & ligand coordinates
                        prot_lig_coords = self.filterCoords(inter, inner_key)                   

                        # create pseudoatoms to show interaction as distance
                        for i, coords in enumerate(prot_lig_coords):
                            cmd.pseudoatom(f"pseudos{i+1}", pos=self.getCoords(coords), state=t+1)

                        # create interactions (as distances)
                        cmd.distance(k, "pseudos1", "pseudos2", state=t+1)
                        cmd.delete("pseudos1, pseudos2")
            
            # create receptor_residues
            # discrete should improve the efficiency
            # cmd.create("receptor_residues", "temp01927824", 1, t+1, discrete=1) TODO -> don't have discrete -> additional sequence is not shown anymore...
            cmd.create("receptor_residues", "temp01927824", 1, t+1)
            

            # clean up
            cmd.delete("temp01927824")
            cmd.delete("tempSeed")

            # create results
            cmd.create("results", f"plip{t}", 1, t+1)
            # set title using filepath (for distance filtering)
            cmd.set_title("results", t+1, str(v["path"]))
            cmd.delete(f"plip{t}")

        # create selection string
        str_unique = " ".join(unique_interactions)

        unique_interactions = set(unique_interactions)
        
        # styling according the different interaction typses
        self.style(unique_interactions)
        
        # group created pymol objects
        cmd.group("interactions", str_unique)
        cmd.group("structures", "results receptor_residues", )
        cmd.group("plip", "interactions structures")
        
        # correct hydrogen
        cmd.remove("hydrogens and plip")
        # cmd.h_add("plip and (not backbone)") 

        # hide labels of lengths
        cmd.hide('labels', 'interactions')
        
        # clean up & set back view
        remaining_interactions = [x for x in all_interactions if x not in unique_interactions]
        [cmd.delete(x) for x in remaining_interactions]
        cmd.set_view(myview)
        cmd.undo_enable()

    def selectAll(self) -> None:
        """
        select / deselect all checkboxes for the interaction types
        """
        # set states according to state of selectAll checkbox
        for k in styles.PYMOL_STYLES.keys():
            self.findChild(QCheckBox, k).setChecked(self.findChild(QCheckBox, "selectAll").isChecked())
        
        # debouncing
        self.debounceTypes.start

    def setInteractionFilter(self) -> None:
        """
        set self.API.interactonFilter based on input
        specify interaction types
        """
        # set interactionFilter
        for k in styles.PYMOL_STYLES.keys():
            self.API.interactionFilter[k] = self.findChild(QCheckBox, k).isChecked()
        
        # functionality of selectAll button (only checked if all are checked)
        if sum(self.API.interactionFilter.values()) == len(self.API.interactionFilter.values()):
            self.findChild(QCheckBox, "selectAll").setChecked(True)
        else:
            self.findChild(QCheckBox, "selectAll").setChecked(False)

        self.showInteractions()

    @staticmethod
    def str2int(value: str) -> int:
        """
        transform an input string into an integer

        exception handling for input residues
        """
        if value.strip():
            try:
                value = int(value)
                if value <= 0:
                    print(f"Warning: Residue {value} selected")
                return value
            except:
                print(f"Incorrect Residue '{value}' selected")

    def setResidueFilter(self, k: str) -> None:
        """
        set self.API.residueFilter based on input
        specify residue for specific interaction

        :param str k: interactiontype e.g., hydrophobic_interactions
        """
        
        # get input
        inputtxt = self.findChild(QLineEdit, f"{k}_res").text()
        
        # drop interactiontypes w/o specified residues from residueFilter
        if inputtxt == "":
            self.API.residueFilter.pop(k, None)
            self.showInteractions()
        
        # add specified residues to residueFilter
        else:
            residues = inputtxt.split(",")
            self.API.residueFilter[k] = [self.str2int(x) for x in residues]
            self.showInteractions()
            

    def setLogicalFilter(self, k: str) -> None:
        """
        set self.API.logicalFilter based on input
        change color of button

        :param str k:   "OR" | "AND"
        """
        self.API.logicalFilter = k
        for key, value in self.logicalPushButtons.items():
            if key == self.API.logicalFilter:
                value.setStyleSheet(styles.LogicalButton)
            else:
                value.setStyleSheet("")

        if k == "OR":
            self.minAndIntCheckBox.setEnabled(False)
            self.minAndIntCheckBox.setChecked(False)
            self.activateSpinBox()
        else:
            self.minAndIntCheckBox.setEnabled(True)

        if self.API.results:
            self.showInteractions()

    def setHideEmpty(self) -> None:
        """
        set self.PI.removeEmpty based on input
        hide states without interactions
        """
        self.API.removeEmpty = self.hideEmptyCheckBox.isChecked()
        self.showInteractions()

    def activateSpinBox(self) -> None:
        """
        activates SpinBox to set minimum number of AND interaction types
        """
        # reset value, if SpinBox is disabled
        if not self.minAndIntCheckBox.isChecked():
            self.minAndInteractions.setValue(0)
        
        # dis/enable SpinBox according to CheckBox
        self.minAndInteractions.setEnabled(self.minAndIntCheckBox.isChecked())

        # update maximum number
        self.minAndInteractions.setMaximum(sum(self.API.interactionFilter.values()))

    def setMinAndInt(self) -> None:
        """
        minimum number of combinations
        """       
        self.API.minAndInteractions = self.minAndInteractions.value()
        self.showInteractions()

    def rightClickedEvent(self) -> None:
        """
        Abd El-Wahed, A. et al. Wasp Venom Biochemical Components and Their Potential in 
        Biological Applications and Nanotechnological Interventions. Toxins 13, 206 (2021).

        General knowledge since the trip to Appenzell.
        """
        # reset GUI, reinitialize PyMOL window and reset self.API
        self.clean_up()

        # return interesting fact & show toxin structure
        print("Mastoparan is a major component of vasp venom. Due to its heat lability,\n" 
              "the symptoms of an insect bite can be minimized by targeted warming.")
        cmd.fetch("2CZP")

    def clearInput(self) -> None:
        """
        reset plugin: clear all inputs  (API is not reset)
        """
        # Distance
        for i in range(2):
            self.findChild(ComboBox, f"Sel_{i+1}").setCurrentText("None")

        self.DistanceLEdit.setText("")
        self.API.distFilter = {}

        # Interactions
        for k in styles.PYMOL_STYLES.keys():
            self.API.interactionFilter[k] = True
            self.findChild(QCheckBox, k).setChecked(True)
            self.findChild(QLineEdit, f"{k}_res").clear()
           
        self.API.residueFilter = {}
        
        self.hideEmptyCheckBox.setChecked(self.API.removeEmpty)
        
        self.setLogicalFilter("OR")

        self.filtflnmLEdit.clear()

        # Show interactions
        if self.API.results:
            self.showInteractions()

    def loadFilter(self) -> None:
        """
        load and apply existing filter.json file
        """
        # reset GUI
        self.clearInput()

        # select json file 
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getOpenFileName(self, "select the filter file", str(self.API.directory / "plip"), 
                                                        "JSON File (*.json, *.JSON)", options=options)
        
        # exception handling if nothing was selected
        if not filepath:
            return

        # load filter in API
        self.API.loadFilter(filepath)
        
        # load filter in GUI
        self.filtflnmLEdit.setText(Path(filepath).stem)

        for k in styles.PYMOL_STYLES.keys():
            self.findChild(QCheckBox, k).setChecked(self.API.interactionFilter[k])

        for k, v in self.API.residueFilter.items():
            self.findChild(QLineEdit, f"{k}_res").insert(", ".join(str(e) for e in v))

        self.setLogicalFilter(self.API.logicalFilter)

        self.hideEmptyCheckBox.setChecked(self.API.removeEmpty)


        if self.API.minAndInteractions:
            self.minAndInteractions.setValue(self.API.minAndInteractions)
            self.minAndInteractions.setEnabled(True)
            self.minAndIntCheckBox.setChecked(True)

        # create Pymol Objects for distance filter
        # Caution: depend on loaded structures, could differ from filter in API 
        if self.API.distFilter:
            self.createPymolSelection()
            self.DistanceLEdit.setText(str(self.API.distFilter["distance"]))

        self.showInteractions()

    def createPymolSelection(self) -> None:
        """
        create pymol objects for the selections of the filter.json file

        called by self.loadFilter()
        """
        # count number of states
        if "results" not in cmd.get_names("all"):
            return

        nstates = cmd.count_states("results")

        # iterate over all states of both selections
        for i in range(2):
            selection = []
            for state in range(1, nstates+1):
                # create selection string for all selected atoms of all states                
                filename = cmd.get_title("results", state)
                if filename in self.API.distFilter[f"Sel_{i+1}"]:
                    for coord in self.API.distFilter[f"Sel_{i+1}"][filename]:
                        selection.append(f"{self.createCoordSelector(coord)} and state {state}")
            strselection = " or ".join(selection)
            # create pymol objects for the selections
            cmd.create(f"Loaded_Selection_{i+1}", strselection)

            # update GUI
            self.updateSelections()

            self.findChild(ComboBox, f"Sel_{i+1}").setCurrentText(f"Loaded_Selection_{i+1}")

    def createCoordSelector(self, coord: tuple) -> str:
        """
        create pymol selection string for an atom by its coordinates

        :param tuple[int] coord:    x, y, z-coordinates for an atom

        >>> createCoordSelector((0, 0, 0))
        "(x > -0.001 and x < 0.001) and (y > -0.001 and y < 0.001) and (z > -0.001 and z < 0.001)"
        """
        # init dim and selector list
        dim = ["x", "y", "z"]
        selector = []

        # iterate over dim and add coordinates
        for d, v in zip(dim, coord):
            selector.append(f"({d} > {self.createLowLim(v)} and {d} < {self.createUpLim(v)})")
        return " and ".join(selector)
    
    @staticmethod
    def createUpLim(v: float, eps: float=0.001) -> float:
        """
        helper function for pymol selection by coordinates
        rounds input value and adds epsilon

        :param float v:     value
        :param float eps:   epsilon which is added to value

        >>> createUpLim(0.250)
        0.251
        """
        return round(v, 3) + eps

    @staticmethod
    def createLowLim(v: float, eps: float=0.001) -> float:
        """
        helper function for pymol selection by coordinates
        rounds input value and subtract epsilon

        :param float v:     value
        :param float eps:   epsilon which is subtracted from value

        >>> createLowLim(0.250)
        0.249
        """
        return round(v, 3) - eps
    
    def saveFilter(self) -> None:
        """
        save filter parameters and results as json file, using filtflnmEdit as file name
        """
        try:
            self.API.saveFilter(self.filtflnmLEdit.text().strip()+".json")
        except:
            print(f"Filter could not be saved as {self.API.plipfolder / self.filtflnmLEdit.text().strip()}.json")

    def showDist(self) -> None:
        """
        Show distance labels for interactions

        called by Menu > Settings > show Distances
        """
        if "interactions" in cmd.get_names("all"):
            if self.MenuBar.actionshowDist.isChecked():
                cmd.show('labels', 'interactions')
            else:
                cmd.hide('labels', 'interactions')

    @staticmethod
    def ngl(selection:str='all') -> None:
        """
        coloring by residue according the NGL Viewer from RCSB Protein Data Bank's website.

        adapted from COLOUR BY RESIDUE - PATRICK BRENNAN, UNIVERSITY OF OXFORD - AUGUST 2020
        https://www.blopig.com/blog/2020/11/pymol-colour-by-residue/
        """
        cmd.select('backbone_carbons','name C+CA')
        cmd.select('ALA','name CB and resn ALA')
        cmd.select('ARG','name CB+CG+CD+CZ and resn ARG')
        cmd.select('ASN','name CB+CG and resn ASN')
        cmd.select('ASP','name CB+CG and resn ASP')
        cmd.select('CYS','name CB and resn CYS')
        cmd.select('GLN','name CB+CG+CD and resn GLN')
        cmd.select('GLU','name CB+CG+CD and resn GLU')
        cmd.select('HIS','name CB+CG+CD2+CE1 and resn HIS')
        cmd.select('ILE','name CB+CG1+CG2+CD1 and resn ILE')
        cmd.select('LEU','name CB+CG+CD1+CD2 and resn LEU')
        cmd.select('LYS','name CB+CG+CD+CE and resn LYS')
        cmd.select('MET','name CB+CG+CE and resn MET')
        cmd.select('PHE','name CB+CG+CD1+CD2+CE1+CE2+CZ and resn PHE')
        cmd.select('PRO','name CB+CG+CD and resn PRO')
        cmd.select('SER','name CB and resn SER')
        cmd.select('THR','name CB+CG2 and resn THR')
        cmd.select('TRP','name CB+CG+CD1+CD2+CE2+CE3+CZ2+CZ3+CH2 and resn TRP')
        cmd.select('TYR','name CB+CG+CD1+CD2+CE1+CE2+CZ and resn TYR')
        cmd.select('VAL','name CB+CG1+CG2 and resn VAL')
        code = {'backbone_carbons':'white','ALA':'lime','ARG':'density','ASN':'deepsalmon','ASP':'warmpink','CYS':'paleyellow','GLN':'tv_red','GLU':'ruby','HIS':'slate','ILE':'forest','LEU':'smudge','LYS':'deepblue','MET':'sand','PHE':'gray40','PRO':'gray20','SER':'tv_orange','THR':'brown','TRP':'palegreen','TYR':'wheat','VAL':'pink'}
        
        cmd.select('none')
        for elem in code:
            cmd.color(code[elem], elem+'&'+selection)
            cmd.delete(elem)

    def clean_up(self) -> None:
        """
        reset GUI, reinitialize PyMOL window and reset self.API
        
        called by Menu > Settings > Clean up
        """
        # reset API
        self.DirectoryLabel.setText("")
        self.API = api.plipAPI(Path(""))
        self.API.pdbonly = self.pdbBtn.isChecked()
        self.MenuBar.submenu_chains.change_chain()
        self.enableEverything(False)
        
        # reset GUI
        self.clearInput()

        # reinitialize PyMOL window
        print("Reinitialize")
        cmd.reinitialize()

# add NGL command to pymol commands
cmd.extend('ngl',InteractionGUI.ngl)
