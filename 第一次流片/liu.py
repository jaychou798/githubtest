
import gdsfactory as gf
import zjupdk
import zjucomponent as zc

# grating = zc.focusin_grating_coupler_y(period=1.033, line_width=0.517, r0=6.1, wgpos=6.1, taper_angle=55.4,
#                                        wg_width=1.334, n_periods=30)

grating=gf.get_component('focusing_grating_dual')

def liu_component(l2,w5):
    c=gf.Component()
    edge = gf.get_component('edge_coupler_overlay', w1=1.134, w2=1.334, w4=0.4, w5=w5, l2=l2)
    edge_r = c << edge
    edge_r1 = c << edge
    grat_r = c << grating
    grat_r1 = c << grating
    grat_r.connect('o1', edge_r.ports['wg_o'])
    edge_r1.connect('fiber_o', edge_r.ports['fiber_o'])
    grat_r1.connect('o1', edge_r1.ports['wg_o'])
    return c

for_output = dict(
    test1 = liu_component(100, 0.4),
    test2 = liu_component(200, 0.4),
    test3 = liu_component(300, 0.4),
    test4 = liu_component(100, 0.35),
    test5 = liu_component(200, 0.35),
    test6 = liu_component(300, 0.35),
)

# if __name__ == "__main__":
#     c = gf.Component()
#
#     ty = 0.0
#     for tc in for_output:
#         c.add_ref(for_output[tc],origin=(0.0,ty))
#         ty = ty +100.0
#
#     c.show()
#     c.write_gds('C:/Users/ROG/Desktop/20230222/liu.gds')
def mycomponent():
    c = gf.Component()
    ty = 0.0
    for tc in for_output:
        c.add_ref(for_output[tc],origin=(0.0,ty))
        ty = ty +100.0
    return c