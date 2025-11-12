import os
import shutil
import numpy as np
import MDAnalysis
from biopandas.pdb import PandasPdb
import warnings
warnings.filterwarnings("ignore")

import time


import json
from pathlib import Path
from collections import defaultdict
from joblib import Parallel, delayed

import xml.etree.cElementTree as ET

from absl import app
from absl import flags
import sys

import subprocess
import concurrent.futures


sys.path.append("/src")

FLAGS = flags.FLAGS
flags.DEFINE_string("dir_path", "/home/pliper/", "")
flags.DEFINE_string("outpath", "plip", "")
flags.DEFINE_boolean('pdbonly', False, 'Process a folder containing only pdb files')
flags.DEFINE_string('intra', '', 'Detect intra-chain interactions')


# TODO: I do not understand the following code fragment / how to find an explanation...

# class FLAGS:
#     dir_path = "/home/pliper/5"
#     outpath = "plip"


def etree_to_dict(t):
    """
    transform ElementTree.Element into defaultdict
    adapted from https://stackoverflow.com/questions/2148119
    # >>> mytree = ET.parse("Path")
    # >>> myroot = mytree.getroot()
    # etree_to_dict(myroot)
    s = '''
            <root>
            <e />
            <e>text</e>
            <e name="value" />
            <e name="value">text</e>
            <e> <a>text</a> <b>text</b> </e>
            <e> <a>text</a> <a>text</a> </e>
            <e> text <a>text</a> </e>
            </root>
        '''
    
    >>> s = '''<root> <e /> <e>text</e> <e name="value" /> <e name="value">text</e> <e> <a>text</a> <b>text</b> </e> <e> <a>text</a> <a>text</a> </e> <e> text <a>text</a> </e> </root>'''
    >>> e = ET.XML(s)
    >>> pprint(etree_to_dict(e))
    {'root': {'e': [None,
                    'text',
                    {'@name': 'value'},
                    {'#text': 'text', '@name': 'value'},
                    {'a': 'text', 'b': 'text'},
                    {'a': ['text', 'text']},
                    {'#text': 'text', 'a': 'text'}]}}
    """
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(("@" + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
              d[t.tag]["#text"] = text
        else:
            d[t.tag] = text
    return d

    
def xmlParse(file):
    """
    Parse xml File and return all data in a defaultdict
    """
    # Create ElementTree (XML Parsing)
    mytree = ET.parse(file)
    myroot = mytree.getroot()

    # Transform ElementTree into defaultdict
    return etree_to_dict(myroot)

def merge_multi_resuts(d):
    merged_dict = {}
    for res in d["report"]["bindingsite"]:
        for key, value in res["interactions"].items():
            if merged_dict.get(key):
                if not value:
                    continue
                for sub_key, sub_value in value.items():
                    merged_dict[key][sub_key] += sub_value
                        
            else:
                if not value:
                    merged_dict[key] = value
                else:
                    for deeper_key, deeper_value in value.items():
                        if type(deeper_value) == list:
                            merged_dict[key] = {deeper_key: deeper_value}
                        else:
                            merged_dict[key] = {deeper_key: [deeper_value]}

    return merged_dict

def edit_distance(string1, string2):
    """
    edit_distance: calculate the minimum number of character insertions, deletions, and substitutions required to change one string to the other
    Ref: https://bit.ly/2Pf4a6Z
    """

    if len(string1) > len(string2):
        difference = len(string1) - len(string2)
        string1 = string1[:-difference]

    elif len(string2) > len(string1):
        difference = len(string2) - len(string1)
        string2 = string2[:-difference]

    else:
        difference = 0

    for i in range(len(string1)):
        if string1[i] != string2[i]:
            difference += 1

    return difference

def write_file(fp):
    file_path = f"{FLAGS.dir_path}/{fp}"
    
    with open(file_path, "r") as r:
        file = r.read()
        
    with open(f"{FLAGS.dir_path}/{FLAGS.outpath}/combined_processed.pdbqt", "a+") as f:
        f.write(file)


## raw input files
def process_template(all_files):
    
    # get pdb file of receptor
    receptor = [x for x in all_files if "_preview.pdb" in x and ".pdbqt" not in x]
    print(receptor)
    
    # get sdf file of ligand
    ligands = [x for x in all_files if "_preview.sdf" in x]

    # combine receptor and ligand
    combined = ""

    for x in receptor:
        combined += FLAGS.dir_path + "/" + x + " "

    for x in ligands:
        combined += FLAGS.dir_path + "/" + x + " "

    # create combined pdb file using openbabel
    command = f"obabel {combined.strip()} -j -d -O {FLAGS.dir_path}/{FLAGS.outpath}/combined.pdb"
    print(command)
    os.system(command)
    # os.system(f"grep 'CONECT' -v {FLAGS.dir_path}/{FLAGS.outpath}/combined.pdb >> {FLAGS.dir_path}/{FLAGS.outpath}/combined_stripped.pdb")
    
    # return name of ligand sdf-file
    return ligands

    
## processed input files
def process_reference(all_files, ligands):
    # get ligand pdbqt files
    processed_ligands = [x.replace(".sdf", ".pdbqt") for x in ligands]
    
    # get flex residue files
    flex_residues = [x for x in all_files if "flex" in x]

    # get pdbqt_vina_out.pdbqt file
    vina_results = [x for x in all_files if "_vina_out.pdbqt" in x][0]

    # edit_distance: calculate the minimum number of character insertions, deletions, and substitutions required to change one string to the other
    multi_ligand_order = np.array([edit_distance(x, vina_results) for x in processed_ligands]).argsort()

    # combine processed files
    try:
        os.remove(f"{FLAGS.dir_path}/{FLAGS.outpath}/combined_processed.pdbqt")
    except:
        pass
    for x in multi_ligand_order:
        write_file(processed_ligands[x])
    for x in flex_residues:
        write_file(x)
    
    # return *_vina_out.pdbqt file
    return vina_results

def get_atom_dict(univ):
    return {x : y.position for (x,y) in enumerate(univ.atoms)}

def get_atom_dict_biopandas(df):
    return {x : y[["x_coord", "y_coord", "z_coord"]].values.astype("float32") for (x,y) in df.iterrows()}

def split_result_states(vina_result_all_states):
    # split / create a pdbqt file for each state

    with open(f"{FLAGS.dir_path}/{vina_result_all_states}", "r") as f:
        file = f.read().split("MODEL")[1:]
    
    processed_files = []
    for t, res in enumerate(file):
        with open(f"{FLAGS.dir_path}/{FLAGS.outpath}/result_{t}.pdbqt", "a+") as handle:
            handle.write("MODEL" + res)
        processed_files.append(f"{FLAGS.dir_path}/{FLAGS.outpath}/result_{t}.pdbqt")
    
    return processed_files


def update_template(result_pdbqt):

    # get file index
    idx = result_pdbqt.split("_")[-1].split(".")[0]
    print(idx)
    
    # transform combined.pdb into pandas dataframe
    pl1 = PandasPdb().read_pdb(f"{FLAGS.dir_path}/{FLAGS.outpath}/combined.pdb")

    # get MDAnalysis...AtomGroup of all atoms w/o hydrogen of combined_processed.pdbqt
    u_pre_dock = MDAnalysis.Universe(f"{FLAGS.dir_path}/{FLAGS.outpath}/combined_processed.pdbqt").select_atoms("not name H*")

    # transform atoms of AtomGroup into dictionary {idx: coordinates}
    pre_dock_atom_dict = get_atom_dict(u_pre_dock)

    # split pdb pandas dataframe into ATOM & HETATM
    file_before_atom_dict = {}
    keys = ["ATOM", "HETATM"]

    # transform atoms into dictionary {idx, coordinates}
    for k in keys:
        file_before_atom_dict[k] = get_atom_dict_biopandas(pl1.df[k])
    
    # create a mapping dictionary {biopandas: MDAnalysis...AtomGroup}
    mapping = {}
    for k in keys:
        mapping_temp = {}
        for k_before, v_before in file_before_atom_dict[k].items():
            for k_pre, v_pre in pre_dock_atom_dict.items():
                if all([np.abs(x-y) <= 0.05 for x,y in zip(v_before, v_pre)]):
                    mapping_temp[k_before] = k_pre
        mapping[k] = mapping_temp

    # get MDAnalysis...AtomGroup of all atoms w/o hydrogen of result_pdbqt
    u_post_dock = MDAnalysis.Universe(f"{result_pdbqt}").select_atoms("not name H*")
    
    # transform atoms of AtomGroup into dictionary {idx: coordinates}
    u_post_dock_atom_dict = get_atom_dict(u_post_dock)

    # add original atom numbering
    for k in keys:
        submap = mapping[k]
        for key_before, k_pre in submap.items():
            pl1.df[k].loc[key_before, ["x_coord", "y_coord", "z_coord"]] = u_post_dock_atom_dict[k_pre]

    # init paths
    outpath = f"{FLAGS.dir_path}/{FLAGS.outpath}/result_{idx}.pdb"
    outpath_h = f"{FLAGS.dir_path}/{FLAGS.outpath}/result_{idx}_h.pdb"
    
    # save result{idx}.pdb with original atom indices
    pl1.to_pdb(path=outpath, 
                records=None, 
                gz=False, 
                append_newline=True)

    # add hydrogens & get resulting file
    os.system(f"obabel {outpath} -O {outpath_h} -h")
    return outpath_h


def getFiles(directory: str) -> list:
    """
    get files from directory
    substitute spaces into underscores of all files of a directory, to prevent problems with openbabel & plip
    
    |- directory
    |----Test Structure 1.pdb
    |----Test_Structure2.pdb

    >>> getFiles(directory)
    ["Test Structure1.pdb", "Test_Structure2.pdb"]

    |- directory
    |----Test_Structure_1.pdb
    |----Test_Structure2.pdb
    """
    filelist = os.listdir(directory)
    
    # exception handling (avoid problems with openbabel and plip)
    for i, filename in enumerate(filelist):
        if " " in filename:
            newfilename = filename.replace(" ", "_")
            filelist[i] = newfilename
            os.rename(f"{FLAGS.dir_path}/{filename}", f"{FLAGS.dir_path}/{newfilename}")
            print(filename, "renamed into", newfilename)
    return filelist

def prepare_plip(pdbonly=False):
    # delete existing output directory
    try:
        shutil.rmtree(f"{FLAGS.dir_path}/{FLAGS.outpath}")
    except:
        pass
    
    # create new output directory
    os.makedirs(f"{FLAGS.dir_path}/{FLAGS.outpath}", exist_ok=True)
    
    # get all input files (all files from directory)
    all_files = getFiles(FLAGS.dir_path)

    if pdbonly:
        # TODO: separate this into a seperate function & parallelize
        # add hydrogens & get resulting file
        all_files_for_plip = []
        for pdbfile in all_files:
            if pdbfile.endswith(".pdb"):
                # init paths
                inpath = f"{FLAGS.dir_path}/{pdbfile}"
                outpath_h = f"{FLAGS.dir_path}/{FLAGS.outpath}/{pdbfile[:-4]}_h.pdb"
                all_files_for_plip.append(outpath_h)
                try:
                    subprocess.run(f"obabel {inpath} -O {outpath_h} -h")
                except:
                    os.system(f"cp {inpath} {outpath_h}")
        return all_files_for_plip

    else:
        # create combined (receptor & ligand) pdb-file & get name of ligand sdf-file
        ligands = process_template(all_files)

        # process and combine files & get *_vina_out.pdbqt 
        vina_result_all_states = process_reference(all_files, ligands)

        # split / create a pdbqt file for each state & get list of all new created files
        result_files = split_result_states(vina_result_all_states)
        
        # TODO: not sure if i got this correct - add hydrogens and keep original atom numbering & get all processed/prepared files
        all_files_for_plip = Parallel(n_jobs=48, backend="multiprocessing")(delayed(update_template)(x) for x in result_files)
        
        return all_files_for_plip


def multi_getresults(f):
    """parallelization of merging different results into one (not sure if needed)"""
    f_short = f.split("/")[-1] #replace(f"{FLAGS.dir_path}/", "")
    result_pdb_name = f_short.split(".")[0]
    
    d = xmlParse(f"{FLAGS.dir_path}/{FLAGS.outpath}/Plip/{result_pdb_name}/report.xml")
    
    # type exception handling
    if type(d["report"]["bindingsite"]) == list:
        interactions = merge_multi_resuts(d)
    else:
        interactions = d["report"]["bindingsite"]["interactions"]
    
    return result_pdb_name, f_short, interactions



def main(argv):
    af = prepare_plip(pdbonly=FLAGS.pdbonly)
    only_one = len(af) == 1
    ff = " ".join(af)       

    #  flag for intra chain interactins
    if FLAGS.intra != "0":
        intra = f"--intra {FLAGS.intra}"
    else:
        intra = ""

    t_start = time.time()
    # send files to plip
    os.system(f"/usr/bin/python3 /src/plip/plipcmd.py -s -f {ff} -o {FLAGS.dir_path}/{FLAGS.outpath}/Plip -x --nohydro {intra}")
    t_end = time.time()

    print(f"Plip took {t_end-t_start:.2f} s")

    result_json = {}

    # get interactions from plip results
    # TODO: use af instead of ff.split(" ")
    filelst = ff.split(" ")

    if only_one:
        futur = [multi_getresults(filelst[0])]
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futur = executor.map(multi_getresults, filelst)

    # merge interactions for resulting json file
    for result_pdb_name, f_short, interactions in futur:
        result_json[result_pdb_name] = {
            "path" : f_short,
            "interactions" : interactions
        }

    # clean up
    print("DELETEING STUFF")
    with open(f"{FLAGS.dir_path}/{FLAGS.outpath}/results.json", "w") as json_file:
        json.dump(result_json, json_file, indent=2)


    shutil.rmtree(f"{FLAGS.dir_path}/{FLAGS.outpath}/Plip")
    
    remove_files = [x for x in os.listdir(f"{FLAGS.dir_path}/{FLAGS.outpath}") if x not in ff]
    remove_files.remove("results.json")
    for file in remove_files:
        os.remove(f"{FLAGS.dir_path}/{FLAGS.outpath}/{file}")



if __name__ == '__main__':
    app.run(main)

