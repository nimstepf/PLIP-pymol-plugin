"""
Application programming interface to access the PLIP features

api.py
"""

import json
import numpy as np
from pathlib import Path
from Bio.PDB import *

import time
import random

try:
    import paramiko
except:
    from pip._internal import main
    main(['install', 'paramiko'])


class plipAPI:
    """
    API to run plip analysis & filter results 
    
    Pymol bug - cannot understand type hinting like the following 
    """
    # interactions: list[str]

    # directory:    Path         # directory with pdb/pdbqt file to analyse (docker result / folder with only pdbfiles...)
    # resultfile:   Path         # results.json
    # results:      dict         # content of resultfile

    # logicalFilter:        Union["OR", "AND"]                          # how to combine interactiontypes & residues
    # interactionFilter:    dict[str, bool]                             # selected interaction types e.g. {'hydrophobic_interactions': True, 'metal_complexes': False, ...}
    # residueFilter:        dict[str, list[int]]                        # selected residues for specific interaction e.g. {"salt_bridges": [165, 211]}
    # removeEmpty:          bool                                        # remove models w/o interactions
    # distFilter:           dict[str, dict[str, list[tuple[float]]]]    # coordinates of selected atoms & distance limit e.g. {"Sel_1": {"xyz1.pdb": [(0, 0, 0), (1, 1, 1), ...], "xyz2.pdb": ...}, ..., "distance": 5.0}
    # biopdb_distFilter     dict[str, tuple[str, int]]                  # identificators for selection (atom / residue)

    # hostname:             str
    # cpu_cnt:              int
    # connect_on_port:      int
    # username:             str
    # password:             str
    

    def __init__(self, directory: Path) -> None:
        
        self.interactions = ['hydrophobic_interactions', 'metal_complexes', 'salt_bridges', 'pi_cation_interactions', 'halogen_bonds', 'pi_stacks', 'water_bridges', 'hydrogen_bonds']

        # initialize directory & results
        self.plipfolder = None
        self.resultfile = None
        self.results = None
        self.intraChain = "0"
        self.pdbonly = False
        self.pdbsave = True
        self.updateDirectory(directory)

        
        # initialize Filter
        self.logicalFilter = "OR"
        self.interactionFilter = {i: True for i in self.interactions}
        self.residueFilter = {}
        self.removeEmpty = True 
        self.minAndInteractions = None
        self.distFilter = {}
        self.biopdb_distFilter = None

        # initialize ssh connection
        self.hostname = "localhost"
        self.cpu_cnt = 4                
        self.connect_on_port = 78       
        self.username = "root"
        self.password = "Sommer2020"

    def updateDirectory(self, newdirectory: Path) -> None:
        """
        update class variables 
            - directory & resultfile
            - results (if file exists)

        :param Path newdirectory: directory with pdb/pdbqt file to analyse (docker result / folder with only pdbfiles...)
        """
        self.directory = newdirectory
        self.plipfolder = newdirectory / "plip"
        self.resultfile = newdirectory / "plip" / "results.json"
        
        # exception handling if result file does not exists
        if self.resultfile.exists():
            self.results = self.loadJSON(self.resultfile)
        else:
            self.results = None

    def run_plip(self) -> None:
        """
        run plip analysis via ssh
        """
        # exception handling if results were already generated
        if self.results:
            print("Done!")

        else:
            print("Running")
            # create plip subfolder for results
            self.plipfolder.mkdir(exist_ok=True)
    
            # collect all relevant files from directory
            files = [file for file in self.directory.glob("*") if file.suffix != ".pse" and file.suffix != ".txt" and file.is_file()]
                       
            # send files to remote and run plip analysis
            self.sync_to_remote_and_task(files, intraChain=self.intraChain, pdbonly=self.pdbonly, pdbsave=self.pdbsave)

            # update results
            if self.resultfile.exists():
                self.results = self.loadJSON(self.resultfile)
                print("Done!")

    def sync_to_remote_and_task(self, files: list, intraChain: str="0", pdbonly: bool=False, pdbsave: bool=True) -> None:
        """
        send files to remote and run plip analysis

        :param list[Path] files: list of files to process
        :param str intraChain    flag to analyse intra chain interactions
        :param bool pdbonly:     flag to differentiate between alphadock results & pdb files only
        :param bool pdbsave:     flag to (not) save processed pdb files
        """
        # initialize SSH client
        try: 
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.load_system_host_keys()
            ssh.connect(hostname=self.hostname, port=self.connect_on_port,
                        username=self.username, password=self.password)
        except:
            print("Error: Could not analyze structures with Plip. Docker Container offline."
                  "Restart docker container on linux with 'docker run -t -d -p 78:22 nimstepf/plip-it:latest' or 'docker restart [container]' or contact Nicolas for help.")
            return

        # create temporary remote folder
        temp_folder = time.ctime().replace(" ", "").replace(":", "") + hex(random.randint(0, 10e8))
        remote_file_path = "/home/pliper/" + temp_folder
        (stdin, stdout, stderr) = ssh.exec_command(command=f"mkdir {remote_file_path}")
        (stdin, stdout, stderr) = ssh.exec_command(f"chmod 0777 {remote_file_path}")
        
        # copy files to remote folder using an ftp client
        ftp_client = ssh.open_sftp()
        for f in files:
            ff_remote = remote_file_path + "/" + f.name
            ftp_client.put(f, ff_remote)
        ftp_client.close()

        # flag to differentiate between alphadock results & pdb files only
        lookup_pdbonly = {True: "--pdbonly", False: ""}
        
        # run plipIt.py script
        (stdin, stdout, stderr) = ssh.exec_command(
            command=f"/root/miniconda3/bin/python /home/pliper/plipIt.py --dir_path {remote_file_path} --intra {intraChain} {lookup_pdbonly[bool(pdbonly or (intraChain!='0'))]}")

        for f in stdout.readlines():
            print(f)

        for f in stderr.readlines():
            print(f)

        # get resulting files
        # create list of files
        (stdin, stdout, stderr) = ssh.exec_command(
            command=f"ls {remote_file_path}/plip")
        retrieve = [x.replace("\n", "") for x in stdout.readlines()]

        # get files via ftp client
        ftp_client = ssh.open_sftp()

        for r in retrieve:
            if pdbsave or r.endswith(".json"):
                fp = remote_file_path + "/plip/" + r
                ff = self.plipfolder / r
                ftp_client.get(fp, ff)

        ftp_client.close()

        # clean up remote folder
        (stdin, stdout, stderr) = ssh.exec_command(
            command=f"rm -rf {remote_file_path}")

    def filterResults(self) -> dict:
        """
        :return:    filtered self.results based on 
                        self.logicalFilter
                        self.interactionFilter
                        self.residueFilter
                        self.removeEmpty
                        self.distFilter        
        """

        filtered_dict = self.results

        # remove results without the selected interaction types
        if self.removeEmpty:
            filtered_dict = self.filterInteractionType(filtered_dict)

        # remove results without the interacting residue
        if self.residueFilter:
            filtered_dict = self.filterInteractionResidue(filtered_dict)

        # distance filter
        if len(self.distFilter) == 3:
            filtered_dict = self.filterDist(filtered_dict)            

        # biopdb distance filter 
        if self.biopdb_distFilter:
            filtered_dict = self.filterBioPDBDistances(filtered_dict)
        
        return filtered_dict

    def filterInteractionType(self, filtered_dict: dict) -> dict:
        """
        reduces the input filtered_dict by removing all items without the selected filter types

        :param dict filtered_dict:      processed self.results dictionary based on applied filters
        :return:                        processed filtered_dict dictionary based on applied filters
        """       
        interactions_to_keep = [x for x in self.interactions if self.interactionFilter[x]]
        
        if self.logicalFilter == "OR":
            return {k: v for k, v in filtered_dict.items() if [True for kt, vt in v["interactions"].items() if kt in interactions_to_keep and vt]}

        elif self.logicalFilter == "AND" and self.minAndInteractions:
            return  {k: v for k, v in filtered_dict.items() if self.checkMin([kt for kt, vt in v["interactions"].items() if vt], interactions_to_keep, self.minAndInteractions)}

        elif self.logicalFilter == "AND":
            return {k: v for k, v in filtered_dict.items() if self.checkAll([kt for kt, vt in v["interactions"].items() if vt], interactions_to_keep)}
            
    @staticmethod
    def checkAll(list1: list, list2: list) -> bool:
        """
        check if the list1 contains all elements of the list2
        adapted from stackoverflow.com/questions/19389490
        """
        return all(item in list1 for item in list2)
    
    @staticmethod
    def checkMin(list1: list, list2: list, minElements: int) -> bool:
        """
        check if the list1 contains a certain number of elements of the list 2
        """
        return sum([item in list1 for item in list2]) >= minElements
    
    def filterInteractionResidue(self, filtered_dict: dict) -> dict:
        """
        reduces the input filtered_dict by removing all items without the specified interacting residue

        :param dict filtered_dict:      processed self.results dictionary based on applied filters
        :return:                        processed filtered_dict dictionary based on applied filters
        """       
        result_to_keep = set()

        for k, v in filtered_dict.items():
            for ks, vs in v["interactions"].items():
                if ks in self.residueFilter and vs:
                    residues = set()
                    
                    # avoid problems with mixed data types
                    iterator = list(vs.values())[0]
                    if isinstance(iterator, dict):
                        iterator = [iterator]

                    # iterate over single interactions from same type
                    for singleint in iterator:
                        # get involved receptor residue
                        residue = int(singleint["resnr"])
                        residues.add(residue)
                       
                    # filter based on residue only
                    if self.logicalFilter == "AND":
                        if all(e in residues for e in self.residueFilter[ks]):
                            result_to_keep.add(k)
                            
                    elif self.logicalFilter == "OR":
                        if [True for x in residues if x in self.residueFilter[ks]]:
                            result_to_keep.add(k)
                
        return {k: v for k, v in filtered_dict.items() if k in result_to_keep}

    def filterDist(self, filter_dict: dict) -> dict:
        """
        reduces the input filter_dict by removing all items breaking distance limit

        :param dict filtered_dict:      processed self.results dictionary based on applied filters
        :return:                        processed filtered_dict dictionary based on applied filters
        """    
        
        files2keep = set()

        # iterate over all atoms from both selections
        for filename, coordinates in self.distFilter[f"Sel_1"].items():
            if filename in self.distFilter[f"Sel_2"]:
                coordinates2 = self.distFilter[f"Sel_2"][filename]

                # check if atoms are within the distance limit
                if self.distCoordFilter(coordinates, coordinates2):
                    files2keep.add(filename)
        
        return {key: value for key, value in filter_dict.items() if value["path"] in files2keep}
        
    def distCoordFilter(self, coords1: list, coords2: list) -> bool:
        """
        iterate over all atom coordinates and return if distance limit is broken

        :param list[tuple[float]] coords1: coordinates of selected atoms (selection 1)
        :param list[tuple[float]] coords2: coordinates of selected atoms (selection 2)
        :return: boolean; True if distance is close enough
        """
        for coord1 in coords1:
            for coord2 in coords2:
                distance = np.linalg.norm(np.array(coord1)-np.array(coord2))
                if distance < self.distFilter["distance"]:
                    return True
        return False

    def saveFilter(self, filename: str) -> None:
        """
        save the different filters in a json file
            - resulting path: self.directroy / plip / filename
            - overwrite: exception handling if file already exists (create filename_{x}.json)

        :param str filename:    filename of the resulting json file e.g., "myfilter.json"
        """
        # merge different filter
        filter = {
            "Interactions": self.interactionFilter,
            "Residues": self.residueFilter,
            "Boolean": self.logicalFilter,
            "Remove": self.removeEmpty,
            "Distance": self.distFilter,
            "Min And": self.minAndInteractions,            
            "Results": list(self.filterResults().keys())
        }

        # save merged filter as json file
        self.saveJSON(filename, filter)
 
    def saveJSON(self, filename: str, data: dict, overwrite=False) -> None:
        """
        save data as json file
            - resulting path: self.directroy / plip / filename
            - overwrite: exception handling if file already exists (create filename_{x}.json)

        :param str filename:    filename of the resulting json file e.g., "myfilter.json"
        :param dict data:       data to save as json {"data": "structured as dictionary"}
        """

        # initialize paths
        filename = Path(filename)
        jsonfile = self.plipfolder / filename
        
        # exception handling if file exists
        i = 1
        while not overwrite:
            if jsonfile.is_file():
                jsonfile = self.plipfolder / (filename.stem + "_" + str(i) + filename.suffix)
                i += 1
            else:
                break

        # save file
        with open(jsonfile, "w") as json_file:
            json.dump(data, json_file, indent=2)

        print(f"File saved as: {jsonfile}")

    def loadJSON(self, filename: str, root=False) -> dict:
        """
        load data from json file (self.plipfolder / filename)

        :param str filename:    filename of the json file e.g., "myfilter.json"
        :param bool root:       exceptionhandling (True for absolute paths, False for relative paths)
        """
        
        # exception handling for relative paths
        if root:
            filepath = Path(filename)
        else:
            filepath = self.plipfolder / filename
        
        # open file
        with open(filepath, "r") as json_file:
            data = json.load(json_file)
        return data

    def loadFilter(self, filename: str) -> None:
        """
        load filter from json file (self.plipfolder / filename)
        
        :param str filename:    filename of the json file e.g., "myfilter.json"
        """
        filter = self.loadJSON(filename)

        self.interactionFilter = filter["Interactions"]
        self.residueFilter = filter["Residues"]
        self.logicalFilter = filter["Boolean"]
        self.removeEmpty = filter["Remove"]
        self.distFilter = filter["Distance"]
        self.minAndInteractions = filter["Min And"]

        print(f"Filter updated from {filename}")

    def filterBioPDBDistances(self, filtered_dict: dict) -> dict:
        """
        reduces the input filter_dict by removing all items breaking distance limit

        :param dict filtered_dict:      processed self.results dictionary based on applied filters
        :return:                        processed filtered_dict dictionary based on applied filters 
        """       
        temp_filtered = {}
        for k, v in filtered_dict.items():
            model = self.getModel(self.plipfolder / v["path"])
            if self.distanceFilter(model):
                temp_filtered[k] = v
        return temp_filtered       

    def distanceFilter(self, model: Model.Model) -> bool:
        """
        parse and transform self.biopdb_distFilter into Atom.Atom/Residue.Residue

        return bool depending on if distance limit is broken
        """
        residues = [self.biopdb_distFilter[f"Sel_{i+1}"] for i in range(2)]
        residues = [self.getResidue(model, residue) if len(residue) == 3 else self.getAtom(model, residue) for residue in residues]
        
        return self.distIterator(residues[0], residues[1], self.biopdb_distFilter["Limit"])

    @staticmethod
    def atom2list(residue):
        """
        helperfunction to process atoms simultaneously with residues

        return list[Atom.Atom] or Residue.Residue
        """
        if isinstance(residue, Atom.Atom):
            return [residue]
        return residue

    def distIterator(self, residue1, residue2, limit: float) -> bool:
        """
        return true if residues/atoms are within the limit 
        
        residue: Bio.PDB.Residue.Residue or Bio.PDB.Atom.Atom
        """
        
        # exception handling to process atoms simultaneously
        residue1 = self.atom2list(residue1)
        residue2 = self.atom2list(residue2)
        
        d = []
        for atom1 in residue1:
            for atom2 in residue2:
                d.append(atom1 - atom2)
        return min(d) < limit

    @staticmethod     
    def getModel(filename: str) -> Model.Model:
        """
        return Bio.PDB.Model.Model for a given pdb file
        """
        parser = PDBParser()
        return parser.get_structure("Structure", filename)[0]

    @staticmethod
    def getAtom(model, atom: tuple) -> Atom.Atom:
        """
        atom: tuple[str, str, int, str, bool]   (tuple with chain, resn, resid, atnm and boolean for atom type [true for hetatm])
        return corresponding atom Bio.PDB.Atom.Atom from model
        """
        chain, resn, resid, atnm, hetatm = atom

        if not hetatm:
            return model[chain][resid][atnm]
        return model[chain][(f"H_{resn}", resid, " ")][atnm]

    @staticmethod
    def getResidue(model, residue: tuple) -> Residue.Residue:
        """
        residue: tuple[str, int, bool]  (tuple with chain, resid, and boolean for atom type [true for hetatm])
        return corresponding residue Bio.PDB.Residue.Residue from model
        """
        chain, resid, hetatm = residue
        
        if hetatm:
            for residue in model[chain]:
                if residue.id[0].startswith("H_") and residue.id[1] == resid:
                    return residue
        return model[chain][resid]

    @staticmethod
    def pdbLine2obj(line: str) -> tuple:
        """
        helperfunction to return the objects to define atom/residue in Bio.PDB
        can be used for getAtom / getResidue

        created using cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        """
        hetatm = line[:6].upper() == "HETATM"
        atnm = line[12:16].strip()
        resn = line[17:20].strip()
        resid = int(line[22:26])
        chain = line[21].strip()
        # Exception handling if no chain is defined
        if not chain: chain = " "

        return chain, resn, resid, atnm, hetatm
