# BetaPLIP

**BetaPLIP** is a Pymol Plugin for the **Protein–Ligand Interaction Profiler (PLIP)** — a tool for automated detection and visualization of non-covalent interactions between proteins and ligands in 3D structures.  
This GUI makes it easier to run PLIP locally and explore results interactively without needing command-line expertise.

## 🧬 Features
- 🖥️ **User-friendly interface** to run PLIP jobs.
- 📁 **Load structures** from files or directories.
- 🔬 **View and analyze** detected protein–ligand interactions.
- 🖼️ **Modern GUI design** with icons and visual feedback.
- 📊 **Integrated result management** for multiple analyses.

## Installation

For the installation and usage of PyMOL, please check out https://www.pymol.org/


You can easily install the GUI by adding the repository to the PyMOL search path:

  1. Click Plugin → Plugin Manager
  2. Go to Settings and click “Add new directory”
  3. Add repository as new directory
  4. Restart PyMOL

The first time you start, a package might download. This should only take a second. Then restart PyMOL again and you should be good to go.


### Docker container to run PLIP analysis locally
Set up the following docker container to run the Plip analysis locally

```docker run -t -d -p 78:22 nimstepf/plip-it:latest```

Don't forget to start the container to analyse structures.


## 📜 License
This project is distributed under the MIT License.
Feel free to modify and adapt it for your own workflow.

## 🧠 Citation
Philipp Schake, Sarah Naomi Bolz, Katja Linnemann, Michael Schroeder, PLIP 2025: introducing protein–protein interactions to the protein–ligand interaction profiler, Nucleic Acids Research, Volume 53, Issue W1, 7 July 2025, Pages W463–W465, https://doi.org/10.1093/nar/gkaf361
Alexander S. Rose, Peter W. Hildebrand, NGL Viewer: a web application for molecular visualization,  Nucleic Acids Research, Volume 43, Issue W1, 1 July 2015, Pages W576–W579. https://doi.org/10.1093/nar/gkv402](https://doi.org/10.1093/nar/gkv402


## ✨ Author
Nicolas Imstepf
🌐 LinkedIn: https://www.linkedin.com/in/nicolas-imstepf/







