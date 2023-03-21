import os.path

from pydantic import BaseModel

import gdsfactory as gf
from gdsfactory.types import Layer, LayerColors, LayerColor, LayerStack, LayerLevel

import sys
sys.path.append(os.path.dirname(__file__))
import zjucomponent as zc

class LayerMapZJU(BaseModel):
    DEVREC: Layer = (68, 0)
    PORT: Layer = (1, 10)  # PinRec optical
    PORTE: Layer = (1, 11)  # PinRec electrical
    SHOW_PORTS: Layer = (1, 12)

    LABEL: Layer = (201, 0)
    LABEL_INSTANCE: Layer = (66, 0)
    TE: Layer = (203, 0)
    TM: Layer = (204, 0)

    WG: Layer = (31, 0)
    WG2_CORE: Layer = (32, 0)
    WG2_CLAD: Layer = (32, 1)
    WG2_REST: Layer = (32, 2)
    METAL1: Layer = (40, 0)
    METAL2: Layer = (41, 0)
    HEATER: Layer = (42, 0)
    HOLE1_CORE: Layer = (50, 0)
    HOLE1_CLAD: Layer = (50, 1)
    HOLE1_REST: Layer = (50, 2)
    HOLE2: Layer = (51, 0)
    HOLE3: Layer = (52, 0)

    class Config:
        frozen = True
        extra = "forbid"
LAYER = LayerMapZJU()

# define layer colors
LAYER_COLORS = LayerColors(layers=dict())

# define cross sections
CROSS_SECTIONS = dict(
    null = gf.partial(gf.cross_section.cross_section, width=0.0, layer=LAYER.PORT),
    strip = gf.partial(gf.cross_section.cross_section, width=1.5, layer=LAYER.WG),
    strip2 = gf.partial(gf.cross_section.cross_section, width=0.28, layer=LAYER.WG2_CORE),
)

# define cells
CELLS_A = dict(
    mmi1x2_y=gf.partial(zc.mmi1x2_y,
                      width=1,
                      width_taper=2,
                      width_mmi=9,
                        length_mmi=20,
                      ),
)

CELLS_name1 = dict(
    mmi1x2_y=gf.partial(zc.mmi1x2_y,
                      width=1.465,
                      width_taper=1.5,
                      width_mmi=9,
                      ),
    mmi1x2_z=gf.partial(zc.mmi1x2_z),
    mmi2x2_y=gf.partial(zc.mmi2x2_y),
    mmi2x2_z=gf.partial(zc.mmi2x2_z),
    eulerbend=gf.partial(zc.eulerbend),

    focusing_grating_air_y=gf.partial(zc.focusin_grating_coupler_y,
                                      wgpos=6.85,
                                      r0=7.15,
                                      taper_angle=55.4,
                                      period=0.995,
                                      line_width=0.398,
                                      wg_width=1.6
                                      ),
    focusing_grating_air_z=gf.partial(zc.focusin_grating_coupler_z,
                                      wgpos=7.466,
                                      r0=7.79,
                                      taper_angle=55.4,
                                      period=0.960,
                                      line_width=0.38,
                                      wg_width=1.6
                                      ),
    focusing_grating_sio2_y=gf.partial(zc.focusin_grating_coupler_y,
                                       wgpos=7.835,
                                       r0=8.18,
                                       taper_angle=55.4,
                                       period=0.965,
                                       line_width=0.382,
                                       wg_width=1.6
                                       ),
    focusing_grating_sio2_z=gf.partial(zc.focusin_grating_coupler_z,
                                       wgpos=10.064,
                                       r0=10.47,
                                       taper_angle=55.4,
                                       period=0.930,
                                       line_width=0.365,
                                       wg_width=1.6
                                       ),
    focusing_grating_dual=gf.partial(zc.focusin_grating_coupler_y,
                                       wgpos=6.1,
                                       r0=6.37,
                                       taper_angle=55.4,
                                       period=1.033,
                                       line_width=0.417,
                                       wg_width=1.334,
                                     n_periods=30
                                       ),

    edge_coupler_overlay = gf.partial(zc.edge_coupler_overlay),
    edge_coupler_double_taper = gf.partial(zc.edge_coupler_double_taper),

    # 双偏振光栅耦合器
    dual_polarization_grating_coupler = gf.partial(
    zc.grating_coupler_y,
    period=1.07,
    width=0.6,
    ),
    # 线性光栅耦合器y传空气包层
    linear_grating_coupler_air_y = gf.partial(
    zc.grating_coupler_y,
    period=0.995,
    width=0.348,
    ),
    # 线性光栅耦合器z传空气包层
    linear_grating_coupler_air_z = gf.partial(
    zc.grating_coupler_z,
    period=0.960,
    width=0.330,
    ),
    # 线性光栅耦合器y传二氧化硅包层
    linear_grating_coupler_sio2_y = gf.partial(
    zc.grating_coupler_y,
    period=0.965,
    width=0.332,
    ),
    # 线性光栅耦合器z传二氧化硅包层
    linear_grating_coupler_sio2_z = gf.partial(
    zc.grating_coupler_z,
    period=0.930,
    width=0.315,
    ),
    bend_dc=gf.partial(zc.bend_dc),
    PSR=gf.partial(zc.PSR),
)

fab_zju=gf.Pdk(
    name="Fab_ZJU",
    cells=CELLS_name1,
    cross_sections=CROSS_SECTIONS,
    layers=LAYER.dict(),
    base_pdk=gf.pdk.GENERIC,
    layer_colors=LAYER_COLORS,
)

fab_A=gf.Pdk(
    name="Fab_A",
    cells=CELLS_A,
    cross_sections=CROSS_SECTIONS,
    layers=LAYER.dict(),
    base_pdk=gf.pdk.GENERIC,
    layer_colors=LAYER_COLORS,
)
fab_zju.activate()


