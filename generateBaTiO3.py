# Author @ Nydia R. Varela-Rosales
# Version v1 2020
# Description: example to visualization of build up structure for quantum mechanical simulation
# Example: visualization of BaTiO3 layers are build-up
# Requires: ase package
import numpy as np
from ase import Atoms, Atom
from ase.io import write, read, vasp
from ase_notebook import AseView, ViewConfig, get_example_atoms
from ase.utils import hsv
from ase.build import surface, add_adsorbate, fcc111
fileIN = "BaTiO.POSCAR"
fileIN1 = "Pt.POSCAR"
positionOnTop = 3.8 # units  armstrongs
vacuumLayer   = 3.8
BaTiO3        = vasp.read_vasp(fileIN)
Pt            = vasp.read_vasp(fileIN1)
slab          = fcc111('Pt', size=(1,1,1),vacuum=vacuumLayer) #Pt
add_adsorbate(slab, BaTiO3, positionOnTop, 'ontop')

config        = ViewConfig()
ase_view      = AseView(config)
s4            = surface(BaTiO3 , (1, 1, 1), 10)
write('slab.png', s4 , rotation='10z,-80x')
write('pt.png', Pt * (3,3,3), rotation='10z,-80x')
write('all.png', slab * (3,3,1), rotation='10z,-80x')

# set up visualization
ase_view = AseView(
    rotations="10z,-80x",
    atom_font_size=16,
    axes_length=30,
    canvas_size=(400, 400),
    zoom=1.2,
    show_bonds=False
)
ase_view.config.uc_dash_pattern           = (.6,.4)
ase_view.config.canvas_color_background   = "blue"
ase_view.config.canvas_background_opacity = 0.2

ase_view.config.atom_lighten_by_depth = 0.7
ase_view.make_svg(slab, center_in_uc=True)

from ase_notebook import concatenate_svgs
from ase_notebook.backend.svg import *
svgs = []
rotationsList = ["45x,45y,45z", "0x", "90x"]
for i in range(len(rotationsList)):
    ase_view.config.rotations = rotationsList[i]
    svgs.append(ase_view.make_svg(slab, center_in_uc=True))
    concatenate_svgs(svgs, max_columns=2, scale=0.5, label=True)
    pdf = svg_to_pdf(ase_view.make_svg(slab, center_in_uc=True), str(i)+"save_file.pdf")

filenameOut = "/media/nydia/sicherung/DFT/BaTiO3/BaTiO3onPt111/"+"POSCAR"
vasp.write_vasp(filenameOut, slab)
