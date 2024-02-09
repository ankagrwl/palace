import os
from math import inf
from pathlib import Path

import gdsfactory as gf
import pyvista as pv
from gdsfactory.components.interdigital_capacitor_enclosed import (
    interdigital_capacitor_enclosed,
)
from gdsfactory.generic_tech import LAYER, get_generic_pdk
from gdsfactory.technology import LayerStack
from gdsfactory.technology.layer_stack import LayerLevel
from IPython.display import display

from gplugins.common.types import RFMaterialSpec
from gplugins.palace import run_capacitive_simulation_palace

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
    )
)
material_spec: RFMaterialSpec = {
    "Si": {"relative_permittivity": 11.45},
    "Nb": {"relative_permittivity": inf},
    "vacuum": {"relative_permittivity": 1},
}



simulation_box = [[-200, -200], [200, 200]]
c = gf.Component("capacitance_palace")
cap = c << interdigital_capacitor_enclosed(
    metal_layer=LAYER.WG, gap_layer=LAYER.DEEPTRENCH, enclosure_box=simulation_box
)
c.add_ports(cap.ports)
substrate = gf.components.bbox(bbox=simulation_box, layer=LAYER.WAFER)
c << substrate
c.plot()
     


# results = run_capacitive_simulation_palace(
#     c,
#     layer_stack=layer_stack,
#     material_spec=material_spec,
#     n_processes=4,
#     element_order=1,
#     simulation_folder=Path(os.getcwd()) / "temporary",
#     mesh_parameters=dict(
#         background_tag="vacuum",
#         background_padding=(0,) * 5 + (700,),
#         port_names=c.ports,
#         default_characteristic_length=200,
#         resolutions={
#             "bw": {
#                 "resolution": 15,
#             },
#             "substrate": {
#                 "resolution": 40,
#             },
#             "vacuum": {
#                 "resolution": 40,
#             },
#             **{
#                 f"bw__{port}": {  # `__` is used as the layer–port delimiter for Palace
#                     "resolution": 20,
#                     "DistMax": 30,
#                     "DistMin": 10,
#                     "SizeMax": 14,
#                     "SizeMin": 3,
#                 }
#                 for port in c.ports
#             },
#         },
#     ),
# )
# display(results)
     
