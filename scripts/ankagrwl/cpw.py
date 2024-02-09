import os
from math import inf
from pathlib import Path

import gdsfactory as gf
import numpy as np
import pyvista as pv
import skrf
from gdsfactory.components.interdigital_capacitor_enclosed import (
    interdigital_capacitor_enclosed,
)
from gdsfactory.generic_tech import LAYER, get_generic_pdk
from gdsfactory.technology import LayerStack
from gdsfactory.technology.layer_stack import LayerLevel
from IPython.display import display
from matplotlib import pyplot as plt

from gplugins.common.types import RFMaterialSpec
from gplugins.palace import run_scattering_simulation_palace

gf.config.rich_output()
PDK = get_generic_pdk()
PDK.activate()

layer_stack = LayerStack(
    layers=dict(
        substrate=LayerLevel(
            layer=LAYER.WAFER,
            thickness=500,
            zmin=0,
            material="Si",
            mesh_order=99,
        ),
        bw=LayerLevel(
            layer=LAYER.WG,
            thickness=200e-3,
            zmin=500,
            material="Nb",
            mesh_order=2,
        ),
        bw_port=LayerLevel(
            layer=LAYER.PORT,
            thickness=200e-3,
            zmin=500,
            material="Nb",
            mesh_order=2,
        ),
    )
)
material_spec: RFMaterialSpec = {
    "Si": {"relative_permittivity": 11.45, "relative_permeability": 1},
    "Nb": {"relative_permittivity": inf, "relative_permeability": 1},
    "vacuum": {"relative_permittivity": 1, "relative_permeability": 1},
}

simulation_box = [[-5000, -5000], [5000, 5000]]
c = gf.Component("scattering_palace")

trace = gf.components.bbox(((-4000, -15), (4000, 15)), layer=LAYER.WG)

c << trace


gnd_1 = gf.components.bbox(((-4000, 33), (4000, 33+800)), layer=LAYER.WG)
gnd_2 = gf.components.bbox(((-4000, -33), (4000, -(33+800))), layer=LAYER.WG)

c << gnd_1

c << gnd_2

c.show()


# Add lumped port rectangles manually, see examples for https://awslabs.github.io/palace/stable/examples/cpw/
lumped_port_1_1 = gf.components.bbox(((-4000, 15), (-4000+18, 15+18)), layer=LAYER.PORT)
lumped_port_1_2 = gf.components.bbox(((-4000, -15), (-4000+18, -(15+18))), layer=LAYER.PORT)
c << lumped_port_1_1
c << lumped_port_1_2
c.add_port("o1_1", lumped_port_1_1.center, layer=LAYER.PORT, width=1, port_type='vertical_te')
c.add_port("o1_2", lumped_port_1_2.center, layer=LAYER.PORT, width=1, port_type='vertical_te')

lumped_port_2_1 = gf.components.bbox(((4000-18, 15), (4000, 15+18)), layer=LAYER.PORT)
lumped_port_2_2 = gf.components.bbox(((4000-18, -15), (4000, -(15+18))), layer=LAYER.PORT)
c << lumped_port_2_1
c << lumped_port_2_2
c.add_port("o2_1", lumped_port_2_1.center, layer=LAYER.PORT, width=1, port_type='vertical_te')
c.add_port("o2_2", lumped_port_2_2.center, layer=LAYER.PORT, width=1, port_type='vertical_te')

substrate = gf.components.bbox(bbox=simulation_box, layer=LAYER.WAFER)
c << substrate
c.show()