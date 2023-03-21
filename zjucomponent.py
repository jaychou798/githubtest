from pydantic import BaseModel
from numpy import *
import gdstk
import gdsfactory as gf
from gdsfactory.component import Component
from gdsfactory.types import Layer, LayerColors, LayerColor, LayerStack, LayerLevel
from gdsfactory.component_layout import _reflect_points
import scipy.integrate as si
from scipy.integrate import quad
import numpy as np
from numpy import ndarray,pi

def offsetx(
        xvals:float
):
    x = [0, 0.23, 0.24, 0.27, 0.33, 0.47, 0.56, 1.40, 4.79, 8.62, 9.62, 12.40, 100]
    y = [0, 0.38, 0.37, 0.50, 0.60, 0.70, 0.80, 1.50, 5.00, 9.00, 10.00, 13.00, 100]
    # 插入点为xvals
    yinterp = interp(xvals, x, y)
    return yinterp

def ellipse_arc(
        a: float,
        b: float,
        theta_min: float,
        theta_max: float,
        angle_step: float = 0.5,
) -> ndarray:
    theta = radians(arange(theta_min, theta_max + angle_step, angle_step))
    xs = a * cos(theta)
    xs = gf.snap.snap_to_grid(xs)
    ys = b * sin(theta)
    ys = gf.snap.snap_to_grid(ys)
    return column_stack([xs, ys])


@gf.cell
def focusin_grating_coupler_y(
        r0: float = 6.6,
        taper_angle: float = 47.66,
        wgpos:float=7.753,
        wg_width: float = 1.2,
        line_width: float = 0.332,
        n_periods: int = 25,
        period: float = 0.965,
)-> Component:
    gap=period-line_width
    r=r0-gap/2
    a=wgpos-r0
    b=wg_width / 2/tan(taper_angle / 2*pi/180)

    tt = [360 * arctan((r - sqrt(
        -a ** 2 * tan(pi * taper_angle / 360) ** 2 - 2 * a * b * tan(pi * taper_angle / 360) ** 2 - b ** 2 * tan(
            pi * taper_angle / 360) ** 2 + r ** 2 * tan(pi * taper_angle / 360) ** 2 + r ** 2)) / (
                                 (a + b - r) * tan(pi * taper_angle / 360))) / pi, 360 * arctan((r + sqrt(
        -a ** 2 * tan(pi * taper_angle / 360) ** 2 - 2 * a * b * tan(pi * taper_angle / 360) ** 2 - b ** 2 * tan(
            pi * taper_angle / 360) ** 2 + r ** 2 * tan(pi * taper_angle / 360) ** 2 + r ** 2)) / ((a + b - r) * tan(
        pi * taper_angle / 360))) / pi]

    rtheta=2 * float(tt[0])

    c = gf.Component()

    taper_arc = ellipse_arc(r, r, -rtheta / 2, rtheta / 2)
    x = -a

    port_position = array((x , 0))
    p0 = port_position + (0, wg_width / 2)
    p1 = port_position + (0, -wg_width / 2)
    pts=vstack([p0, p1, taper_arc])

    c.add_polygon(pts, layer=gf.get_layer('WG'))
    #c << gf.components.rectangle((42.0,42.0), layer=gf.get_layer('WG2'), centered=True).move((21.0+a,0))

    for p in range(n_periods):
        bend = gf.components.bend_circular(radius=r+gap+line_width/2 + p * period, angle=rtheta,
                                           width=line_width, cross_section='cross_section')
        ref = c << bend
        ref.rotate(angle=90 - rtheta / 2, center=(0, r+gap+line_width/2 + p * period))
        ref.movey(-(r+gap+line_width/2 + p * period))
    c.add_port('o1', center=(-a, 0), orientation=180, width=wg_width,
               layer=gf.get_layer('PORT'))

    c.absorb(ref)

    return c

@gf.cell
def focusin_grating_coupler_z(**kwargs):
    c=gf.Component()
    ref=c<<focusin_grating_coupler_y(**kwargs)
    ref.rotate(-90)
    c.add_port('o1',port=ref.ports['o1'])
    return c

@gf.cell
def edge_coupler_overlay_discrepent(
    w1: float = 1.0,
    w2: float = 1.2,
    w3: float = 5.0,
    w4: float = 0.28,
    l1: float = 70.0,
    l2: float = 70.0,
    dis_y: float = 0.0,
) -> Component:

    c = gf.Component()

    s=array([(0.0,-w2/2.0),(l2,-w1/2.0),(l2+l1*2,-w1/2.0)])
    s=append(s,_reflect_points(flip(s,0),(0.0,0.0),(1.0,0.0)),0)
    c.add_polygon(s,layer=gf.get_layer('WG'))

    s=array([(0.0,-w3/2.0),(l2-w4*l1/w1,-w1/2.0-w4),(l2,-w1/2.0),(l2+l1,w1/2.0),(l2+l1*2,w1/2.0+w1),(l2+l1*2,w1/2.0+w1+w4),
             (l2+l1,w1/2.0+w4),(l2,w1/2.0+w4*2),(0.0,w3/2.0)])
    c.add_polygon(s+(0,dis_y), layer=gf.get_layer('WG2'))
    P = gf.Path()
    P.append(gf.path.arc(radius=30, angle=20))
    c << gf.path.extrude(P.move(destination=(l2+l1*2,w1/2.0+w1+w4/2.0+dis_y)), layer=gf.get_layer('WG2'), width=w4)

    c.add_port('wg_o',(0,0),orientation=180.0,width=w2,layer=gf.get_layer('WG'))
    c.add_port('fiber_o', (l2+l1*2, 0),orientation=0.0,width=w1,layer=gf.get_layer('WG2'))

    return c

@gf.cell
def edge_coupler_overlay(
    w1: float = 1.0,
    w2: float = 1.2,
    w3: float = 5.0,
    w4: float = 0.28,
    w5: float = 0.28,
    l1: float = 80.0,
    l2: float = 200.0,
    l3: float = 250.0,
    dis_y: float = 0.0,
    with_second_taper: bool = True,
) -> Component:

    c = gf.Component()

    s=array([(0.0,w3/2.0),(l1,w1/2.0),(l1+l2,w1/2.0)])
    s=append(s,_reflect_points(flip(s,0),(0.0,0.0),(1.0,0.0)),0)
    c.add_polygon(s,layer=gf.get_layer('WG2_CORE'))

    s=array([(0.0,-w2/2.0),(l1,-w4/2.0),(l1+l2,w1),(l1+l2+20.0,w1*1.5),(l1+l2+20.0,w1*1.5+w4),(l1+l2,w1+w4),
             (l1,w4/2.0),(0.0,w2/2.0)])
    c.add_polygon(s+(0,dis_y), layer=gf.get_layer('WG'))
    P = gf.Path()
    P.append(gf.path.arc(radius=30, angle=20))
    c << gf.path.extrude(P.move(destination=(l1+l2+20.0,w1*1.5+w4*0.5+dis_y)), layer=gf.get_layer('WG'), width=w4)

    c.add_port('wg_o',(0,dis_y),orientation=180.0,width=w2,layer=gf.get_layer('WG'))
    c.add_port('fiber_o', (l1+l2, 0),orientation=0.0,width=w1,layer=gf.get_layer('WG2_CORE'))

    if(with_second_taper):
        tc = gf.components.taper(l3,w1,w5,with_bbox=False,cross_section='strip2')
        tc_r = c << tc
        tc_r.connect('o1',c.ports['fiber_o'])
        c.absorb(tc_r)
        c.ports['fiber_o'].move((l3+1.0,0.0))

    c.add_ref(gf.components.rectangle((30.0,5.0), layer=gf.get_layer('HOLE2')),columns=int((l2+l3)/40.0),rows=1,
              spacing=(-40.0,0.0),origin=(l1+l2+l3-41.0,-10.0))
    c.add_ref(gf.components.rectangle((30.0,5.0), layer=gf.get_layer('HOLE3')),columns=int((l2+l3)/40.0),rows=1,
              spacing=(-40.0,0.0),origin=(l1+l2+l3-41.0,-10.0))

    tp=gf.Port('',0.0,(l1,0.0),cross_section=gf.get_cross_section('null'))
    c.add_ref(gf.components.rectangle((l2+l3+1,5.0),layer=gf.get_layer('HOLE1_CORE'))).connect('e3',tp)
    c.add_ref(gf.components.rectangle((l2+l3+1,20.0),layer=gf.get_layer('HOLE1_CLAD'))).connect('e3',tp)

    tp = gf.Port('', 0.0, (0.0, 0.0), cross_section=gf.get_cross_section('null'))
    c.add_ref(gf.components.rectangle((l1+l2+l3+1,125.0),layer=gf.get_layer('WG2_CLAD'))).connect('e3',tp)

    return c

@gf.cell
def edge_coupler_double_taper(
    w1: float = 1.2,
    w2: float = 5.0,
    w3: float = 0.25,
    w4: float = 1.5,
    w6: float = 1.0,
    l1: float = 100.0,
    l3: float = 20.0,
    dis_y: float = 0.0,
) -> Component:

    c = gf.Component()

    s=array([(0.0,-w1/2.0),(l1,-w3/2.0),(l1,w3/2.0),(0.0,w1/2.0)])
    c.add_polygon(s,layer=gf.get_layer('WG'))

    s=array([(0.0,-w2/2.0),(l1,-w4/2.0),(l1+l3,-w6/2.0)])
    s=append(s,_reflect_points(flip(s,0),(0.0,0.0),(1.0,0.0)),0)
    c.add_polygon(s+(0,dis_y), layer=gf.get_layer('WG2'))

    c.add_port('wg_o',(0.0,0.0),orientation=180.0,width=w1,layer=gf.get_layer('WG'))
    c.add_port('fiber_o', (l1+l3,dis_y),orientation=0.0,width=w6,layer=gf.get_layer('WG2'))

    return c
@gf.cell
def bend_dc(
        r: float = 140,
        width: float = 0.8,
        gap: float = 0.6,
        theta1: float = 10,
        theta2: float = 45
):
    c = gf.Component()
    r2 = r + gap+width
    c1 = c<<gf.components.bend_circular(radius=r2, angle=theta1, width=width)#下段圆弧
    c1.rotate(angle=-theta1/2, center=(0,r2))
    c2 = c<<gf.components.bend_circular(radius=r, angle=theta2, width=width)#上段圆弧
    c2.rotate(angle=-theta2/2,center=(0,r))
    c3 = c<<gf.components.bend_circular(radius=r2 , angle=theta1/2, width=width)#下段圆弧衔接
    c5 = c << gf.components.bend_circular(radius=r2 , angle=theta1 / 2, width=width)
    c4=c<<gf.components.bend_circular(radius=r*3, angle=theta2/2, width=width)#上段圆弧衔接
    c6 = c << gf.components.bend_circular(radius=r * 3, angle=theta2 / 2, width=width)

    c3.mirror(p1=[0, 0], p2=[1, 0])
    c3.connect('o1', c1.ports['o2'])
    c5.mirror(p1=[0, 0], p2=[1, 0])
    c5.connect('o2', c1.ports['o1'])

    c2.move(origin=[0, 0], destination=[0, width + gap])

    c4.mirror(p1=[0, 0], p2=[1, 0])
    c4.connect('o1', c2.ports['o2'])
    c6.mirror(p1=[0, 0], p2=[1, 0])
    c6.connect('o2', c2.ports['o1'])

    stai1 = c << gf.components.straight(length=(c2.xsize * 4 - c1.xsize * 5 / 4)/2, width=width)
    stai1.connect('o1', c3.ports['o2'])
    stai2 = c << gf.components.straight(length=(c2.xsize * 4 - c1.xsize * 5 / 4)/2, width=width)
    stai2.connect('o2', c5.ports['o1'])

    c.add_port('o1', port=stai2.ports['o1'], layer=gf.get_layer('PORT'))
    c.add_port('o2', port=stai1.ports['o2'], layer=gf.get_layer('PORT'))
    c.add_port('o3', port=c6.ports['o1'], layer=gf.get_layer('PORT'))
    c.add_port('o4', port=c4.ports['o2'], layer=gf.get_layer('PORT'))
    c.absorb(c1)
    c.absorb(c2)
    c.absorb(c3)
    c.absorb(c4)
    c.absorb(c5)
    c.absorb(c6)
    return c

@gf.cell
def PSR(
        width1:float=1.2,
        width2:float=2.1,
        width3:float=2.6,
        l1:float=800,
        l2:float=50,
        width4: float = 0.5,
        width5: float = 0.9,
        width6: float = 2.1,
        gap1: float = 1.5,
        gap2: float = 0.6,
        gap3:float =0.6,
        gap4:float=1.5,
        l3: float = 100,
        l4: float = 200,
        l5: float = 100,
        width7: float = 2.4,
        width8: float = 5.92,
        l6: float = 20,
        l7: float = 53,
        **kwargs
):

    #psr1
    psr1 = gf.Component()
    # c1 = c << gf.components.taper(width1=0.9, width2=1.2, length=50)
    psr12 = psr1 << gf.components.taper(width1=width1, width2=width2, length=l1)
    psr13 = psr1 << gf.components.taper(width1=width2, width2=width3, length=l2)
    # c2.connect('o1', c1.ports['o2'])
    psr13.connect('o1', psr12.ports['o2'])
    psr1.add_port('o1', port=psr12.ports['o1'])
    psr1.add_port('o2', port=psr13.ports['o2'])

    #psr2
    psr2 = gf.Component()
    # psr21 = psr2 << gf.components.taper(width1=width3, width2=width6, length=200)  # 中间下段
    poly=gf.Component()
    poly.add_polygon([(0,width3/2),(l4,width3/2),(l4,width3/2-width6),(0,-width3/2)],layer=gf.get_layer('WG'))
    poly.add_port(name='o1', center=[0,0], width=width3, orientation=180, layer=gf.get_layer('WG'))
    poly.add_port(name='o2', center=[l4,width3/2-width6/2], width=width6, orientation=0, layer=gf.get_layer('WG'))
    psr21=psr2<<poly
    poly1 = gf.Component()
    poly1.add_polygon(
        [(-l3, width3/2 + gap1), (0, width3/2+ gap2), (l4, width3/2 + gap3), (l4 + l5, width3/2 + gap4), (l4 + l5, width3/2 + gap4 + width5),
         (l4, width3/2 + gap3 + width5), (0, width3/2 + gap2 + width4), (-l3, width3/2 + gap1 + width4)], layer=gf.get_layer('WG'))  # 中间上段
    poly1.add_port(name='o1', center=[-l3, width3/2 + gap4 + width4 / 2], width=width4, orientation=180, layer=gf.get_layer('PORT'))
    poly1.add_port(name='o2', center=[l4 + l5, width3/2 + gap4 + width5 / 2], width=width5, orientation=0, layer=gf.get_layer('PORT'))
    p1 = psr2 << poly1
    psr22 = psr2 << gf.components.straight(width=width3, length=l3)  # 左边下段
    psr22.connect('o2', psr21.ports['o1'])
    psr23 = psr2 << gf.components.straight(width=width6, length=l5)  # 右边下段
    psr23.connect('o1', psr21.ports['o2'])
    bend = gf.components.bend_circular(angle=40, radius=90, width=width5)
    ref1 = psr2 << bend
    ref1.connect('o1', poly1.ports['o2'])
    ref2 = psr2 << bend
    ref2.mirror(p1=(0, 0), p2=(1, 0))
    ref2.connect('o1', ref1.ports['o2'])

    psr2.add_port('o1', port=psr22.ports['o1'])
    psr2.add_port('o2', port=psr23.ports['o2'])
    psr2.add_port('o3', port=ref2.ports['o2'])

    #psr3
    psr3 = gf.Component()
    psr31 = psr3 << gf.components.taper(width1=width6, width2=width7, length=l6)
    psr32 = psr3 << gf.components.straight(length=l7, width=width8)
    psr32.connect('o1', psr31.ports['o2'])
    psr33 = psr3 << gf.components.straight(length=10, width=width7)
    psr33.connect('o1', psr32.ports['o2'])

    psr3.add_port('o1', port=psr31.ports['o1'])
    psr3.add_port('o2', port=psr33.ports['o2'])

    c = gf.Component()
    c1 = c << psr1
    c2 = c << psr2
    c3 = c << psr3
    c2.connect('o1', c1.ports['o2'])
    c3.connect('o1', c2.ports['o2'])
    c.add_port('o1', port=c1.ports['o1'])
    c.add_port('o2', port=c3.ports['o2'])
    c.add_port('o3', port=c2.ports['o3'])
    c.absorb(c1)
    c.absorb(c2)
    c.absorb(c3)

    return c

@gf.cell
def grating_coupler_y(
        n: int = 25,#周期数
        wg_width: float = 13,#波导宽度
        period: float = 0.954,#周期间隔
        width: float = 0.381,#光栅齿宽
):
    c = gf.Component()
    polygon = gf.components.rectangle(size=(width, wg_width), layer=gf.get_layer('WG'))
    c.add_array(polygon, columns=n, rows=1, spacing=[period,0])
    gap=period-width
    polygon1 = gf.components.rectangle(size=(4, wg_width), layer=gf.get_layer('WG'))
    polygon2 = gf.components.rectangle(size=(10, wg_width), layer=gf.get_layer('WG'))
    p1=c<<polygon1
    p1.move(origin=p1.center,destination=[-gap-2,wg_width/2])
    p2=c<<polygon2
    p2.move(origin=p2.center,destination=[n*period+5,wg_width/2])
    c.add_port('o1',center=[n*period+10,wg_width/2],orientation=0,width=wg_width,layer=gf.get_layer('PORT'))
    c.absorb(p1)
    c.absorb(p2)
    # ref=t<<c
    # t.absorb(ref)
    # t.add_port('o1', center=[n * period + 10, wg_width / 2], orientation=0, width=wg_width, layer=gf.get_layer('WG'))

    return c

@gf.cell
def grating_coupler_z(
        period:float=0.923,
        width:float=0.369,
        **kwargs
):
    c=gf.Component()
    ref=c<<grating_coupler_y(period=period,width=width,**kwargs)
    ref.rotate(-90)
    c.add_port('o1',port=ref.ports['o1'])
    c.absorb(ref)
    return c

@gf.cell
def eulerbend(
        R_max:float=200,
        R_min:float=37,
        width:float=1.2,
):
    l0=np.pi/2/(1/R_max+1/R_min)
    a=(l0/(1/R_min-1/R_max))**0.5
    num=50
    x=np.zeros(num)
    y = np.zeros(num)
    lmax=l0
    ls=np.linspace(0,lmax,num)
    for i in range(num):
        l=ls[i]
        y[i]=a*si.quad(lambda t:np.sin(t**2/2+a*t/R_max),0,l/a)[0]
        x[i] = a * si.quad(lambda t: np.cos(t ** 2 / 2 + a * t / R_max), 0, l / a)[0]
    x1 = np.zeros(99)
    x1[0:50] = x
    x1[50:99] = -y[-1:0:-1] + x[-1] + y[-1]
    y1 = np.zeros(99)
    y1[0:50] = y
    y1[50:99] = -x[-1:0:-1] + x[-1] + y[-1]
    points = np.column_stack([x1, y1])
    # c = gf.Component()
    p=gf.Path(points)
    c=gf.path.extrude(p,layer=gf.get_layer('WG'),width=width)
    return c

@gf.cell
def mmi1x2_y(
    width: float = 1.5,
    width_taper: float = 3,
    length_taper: float = 20.0,
    length_mmi: float = 57,
    width_mmi: float = 9,
    gap_mmi: float = 1.8,
):
    c = Component()
    gap_mmi = gf.snap.snap_to_grid(gap_mmi, nm=2)
    w_mmi = width_mmi
    w_taper = width_taper

    taper = gf.get_component(
        gf.components.taper,
        length=length_taper,
        width1=width,
        width2=w_taper,

    )
    a = gap_mmi / 2 + width_taper / 2
    mmi = c << gf.get_component(
        gf.components.straight, length=length_mmi, width=w_mmi
    )

    ports = [
        gf.Port(
            "o1",
            orientation=180,
            center=(0, 0),
            width=w_taper,
            layer=gf.get_layer('PORT'),

        ),
        gf.Port(
            "o2",
            orientation=0,
            center=(+length_mmi, +a),
            width=w_taper,
            layer=gf.get_layer('PORT'),

        ),
        gf.Port(
            "o3",
            orientation=0,
            center=(+length_mmi, -a),
            width=w_taper,
            layer=gf.get_layer('PORT'),

        ),
    ]

    for port in ports:
        taper_ref = c << taper
        taper_ref.connect(port="o2", destination=port)
        c.add_port(name=port.name, port=taper_ref.ports["o1"])
        c.absorb(taper_ref)

    c.absorb(mmi)
    return c

@gf.cell
def mmi1x2_z(
    **kwargs
):
    c=gf.Component()
    ref=c<<mmi1x2_y(length_mmi=52,**kwargs)
    ref.rotate(-90)
    c.add_port('o1',port=ref.ports['o1'])
    c.add_port('o2', port=ref.ports['o2'])
    c.add_port('o3', port=ref.ports['o3'])
    return c

@gf.cell
def mmi2x2_y(
    width: float = 1.5,
    width_taper: float = 1.5,
    length_taper: float = 10.0,
    length_mmi: float = 133,
    width_mmi: float = 12,
    gap_mmi: float = 2.5,
):
    c = gf.Component()
    gap_mmi = gf.snap.snap_to_grid(gap_mmi, nm=2)
    w_mmi = width_mmi
    w_taper = width_taper

    taper = gf.get_component(
        gf.components.taper,
        length=length_taper,
        width1=width,
        width2=w_taper,

    )

    a = gap_mmi / 2 + width_taper / 2
    mmi = c << gf.get_component(
        gf.components.straight, length=length_mmi, width=w_mmi
    )

    ports = [
        gf.Port("o1", orientation=180, center=(0, -a), width=w_taper, layer=gf.get_layer('PORT'),),
        gf.Port("o2", orientation=180, center=(0, +a), width=w_taper, layer=gf.get_layer('PORT'),),
        gf.Port(
            "o3",
            orientation=0,
            center=(length_mmi, +a),
            width=w_taper,
            layer=gf.get_layer('PORT'),

        ),
        gf.Port(
            "o4",
            orientation=0,
            center=(length_mmi, -a),
            width=w_taper,
            layer=gf.get_layer('PORT'),

        ),
    ]

    for port in ports:
        taper_ref = c << taper
        taper_ref.connect(port="o2", destination=port)
        c.add_port(name=port.name, port=taper_ref.ports["o1"])
        c.absorb(taper_ref)
    c.absorb(mmi)
    return c

@gf.cell
def mmi2x2_z(
        **kwargs
):
    c = gf.Component()
    ref = c << mmi2x2_y( **kwargs,length_mmi=120,gap_mmi=2.9)
    ref.rotate(-90)
    c.add_port('o1', port=ref.ports['o1'])
    c.add_port('o2', port=ref.ports['o2'])
    c.add_port('o3', port=ref.ports['o3'])
    c.add_port('o4', port=ref.ports['o4'])
    return c
@gf.cell
def wg_cross(
        wg:float=1,
        width:float=3.64,
        l_mmi:float=30,
        l_taper:float=10,
        offset:bool=False
):
    if offset:
        width=offsetx(width)
        wg=offsetx(wg)
    c=gf.Component()
    mmi=c<<gf.components.taper(width1=width,width2=width,length=l_mmi)
    taper1=c<<gf.components.taper(width1=wg,width2=width,length=l_taper)
    taper2 = c << gf.components.taper(width1=wg, width2=width, length=l_taper)
    mmi.center=(0,0)
    taper1.connect('o2',mmi.ports['o1'])
    taper2.connect('o2', mmi.ports['o2'])
    c.add_port('o1',port=taper1.ports['o1'])
    c.add_port('o2', port=taper2.ports['o1'])
    cross=gf.Component()
    ref1=cross<<c
    ref2=cross<<c
    ref1.rotate(45)
    ref2.rotate(-45)
    cross.add_port('o1',port=ref1.ports['o1'])
    cross.add_port('o2', port=ref1.ports['o2'])
    cross.add_port('o3', port=ref2.ports['o1'])
    cross.add_port('o4', port=ref2.ports['o2'])
    cross.absorb(ref1)
    cross.absorb(ref2)
    return cross


def pack(
        com: Component,
        couplery: Component=grating_coupler_y(),
        couplerz: Component=grating_coupler_z(),
        index:float=250/11.5
):
    c = gf.Component()
    ref = c << com
    for i in range(len(com.ports)):

        if ref.ports['o{}'.format((i + 1))].orientation > -5 and ref.ports['o{}'.format((i + 1))].orientation< 5:
            refcoupler = c << couplery
            reftaper = c << gf.components.taper(width1=ref.ports['o{}'.format((i + 1))].width, width2=couplery.ports['o1'].width,
                                                length=10+abs(ref.ports['o{}'.format((i + 1))].width-couplery.ports['o1'].width)*index)
            reftaper.connect('o1', ref.ports['o{}'.format((i + 1))])
            refcoupler.connect('o1', reftaper.ports['o2'])
            c.absorb(refcoupler)
            c.absorb(reftaper)
        if ref.ports['o{}'.format((i + 1))].orientation > 85 and ref.ports['o{}'.format((i + 1))].orientation< 95:
            refcoupler = c << couplerz
            reftaper = c << gf.components.taper(width1=ref.ports['o{}'.format((i + 1))].width, width2=couplerz.ports['o1'].width,
                                                length=10+abs(ref.ports['o{}'.format((i + 1))].width-couplerz.ports['o1'].width)*index)
            reftaper.connect('o1', ref.ports['o{}'.format((i + 1))])
            refcoupler.connect('o1', reftaper.ports['o2'])
            c.absorb(refcoupler)
            c.absorb(reftaper)
        if ref.ports['o{}'.format((i + 1))].orientation > 175 and ref.ports['o{}'.format((i + 1))].orientation < 185:
            refcoupler = c << couplery
            reftaper = c << gf.components.taper(width2=ref.ports['o{}'.format((i + 1))].width, width1=couplery.ports['o1'].width,
                                                length=10+abs(ref.ports['o{}'.format((i + 1))].width-couplery.ports['o1'].width)*index)
            reftaper.connect('o2', ref.ports['o{}'.format((i + 1))])
            refcoupler.connect('o1', reftaper.ports['o1'])
            c.absorb(refcoupler)
            c.absorb(reftaper)
        if ref.ports['o{}'.format((i + 1))].orientation > 265 and ref.ports['o{}'.format((i + 1))].orientation < 275:
            refcoupler = c << couplerz
            reftaper = c << gf.components.taper(width1=ref.ports['o{}'.format((i + 1))].width, width2=couplerz.ports['o1'].width,
                                                length=10+abs(ref.ports['o{}'.format((i + 1))].width-couplerz.ports['o1'].width)*index)
            reftaper.connect('o1', ref.ports['o{}'.format((i + 1))])
            refcoupler.connect('o1', reftaper.ports['o2'])
            c.absorb(refcoupler)
            c.absorb(reftaper)

    c.absorb(ref)
    return c

def dpi():
    d = 0.005
    width = 0.3
    gap = 0.3
    c = gf.Component()
    polygon = gf.components.rectangle(size=(0.4, 1), layer=gf.get_layer('WG'))
    polygon1 = gf.components.rectangle(size=(0.4, 3), layer=gf.get_layer('WG'))


    for i in range(130):
        p = c << gf.components.rectangle(size=(0.3 + i * d, 13), layer=gf.get_layer('WG'))
        p.movex((0.6 + i * d) * i)
    for i in range(150):
        if i % 5 == 0:
            ref = c << polygon1
        else:
            ref = c << polygon
        ref.movex(1 * i)
        ref.movey(14)

    return c


def wg_offset():
    c = gf.Component()
    for i in range(30):
        p = c << gf.components.rectangle(size=(200, 0.25 + i * 0.05), layer=gf.get_layer('WG'))
        p.movey(8 * i)
    return c


def wg_cross(
        wg:float=1,
        width:float=3.64,
        l_mmi:float=30,
        l_taper:float=10,
        offset:bool=False
):
    if offset:
        width=offsetx(width)
        wg=offsetx(wg)
    c=gf.Component()
    mmi=c<<gf.components.taper(width1=width,width2=width,length=l_mmi)
    taper1=c<<gf.components.taper(width1=wg,width2=width,length=l_taper)
    taper2 = c << gf.components.taper(width1=wg, width2=width, length=l_taper)
    mmi.center=(0,0)
    taper1.connect('o2',mmi.ports['o1'])
    taper2.connect('o2', mmi.ports['o2'])
    c.add_port('o1',port=taper1.ports['o1'])
    c.add_port('o2', port=taper2.ports['o1'])
    cross=gf.Component()
    ref1=cross<<c
    ref2=cross<<c
    ref1.rotate(45)
    ref2.rotate(-45)
    cross.add_port('o1',port=ref1.ports['o1'])
    cross.add_port('o2', port=ref1.ports['o2'])
    cross.add_port('o3', port=ref2.ports['o1'])
    cross.add_port('o4', port=ref2.ports['o2'])
    cross.absorb(ref1)
    cross.absorb(ref2)
    return cross

def rightbend(
        width:float=1,
):
    c = gf.Component()
    rightpoints = np.loadtxt('points.txt')
    r = gf.Path(rightpoints)
    you = c << gf.path.extrude(r, layer=gf.get_layer('WG'), width=width)
    c.add_port('o1',port=you.ports['o2'])
    c.add_port('o2', port=you.ports['o1'])
    return c
def spiral(
        width:float=1,
        n_loop:int=20,
        mid_length:float=100,
):
    c = gf.Component()
    w = 1.5 + 5
    length=0
    s1 = c << gf.components.straight(width=width,length=mid_length)
    s2 = c << gf.components.straight(width=width,length=mid_length)
    s0 = c << gf.components.straight(width=width,length=mid_length)
    for i in range(1, n_loop + 1):
        if i == 2:
            Rmin = 37 * (60 + w * 0.5) / 59.693
            Rmax = 200 * (60 + w * 0.5) / 59.693
        else:
            Rmin = 37 * (60 - 3 * w / 2 + w * i) / 59.693
            Rmax = 200 * (60 - 3 * w / 2 + w * i) / 59.693

        L1 = 0
        L2 = 0
        L0 = np.pi / 2 / (1 / Rmin + 1 / Rmax)
        A = np.sqrt(L0 / (1 / Rmin - 1 / Rmax))
        num = 300 + 10 * i
        x = np.zeros(num)
        y = np.zeros(num)
        Lmax = L0

        def fun_y(x):
            return np.sin(x ** 2 / 2 + A * x / Rmax)

        def fun_x(x):
            return np.cos(x ** 2 / 2 + A * x / Rmax)

        Ls = np.linspace(0, Lmax, num=num)
        for ii in range(num):
            L = Ls[ii]
            x[ii] = A * quad(fun_x, 0, L / A)[0]
            y[ii] = A * quad(fun_y, 0, L / A)[0]

        x0 = np.concatenate((x, -y[-2::-1] + x[-1] + y[-1]))
        y0 = np.concatenate((y, -x[-2::-1] + x[-1] + y[-1]))

        x1 = np.flip(x0)
        y1 = 2 * (L1 / 2.0 + np.max(y0)) - np.flip(y0)
        aa = np.vstack((np.hstack((x0, x1)), np.hstack((y0, y1)))).T

        aa[:, 1] = aa[:, 1] - y0[-1]
        rightpoints = aa
        r = gf.Path(rightpoints)
        aa[:, 0] = -aa[:, 0]
        leftpoints = aa
        l = gf.Path(leftpoints)
        if i == 1:
            you = c << rightbend(width=width)
            zuo = c << rightbend(width=width)
            zuo.mirror(p1=(0, 0), p2=(0, 1))
            zuo.connect('o2', s0.ports['o1'])
            s1.connect('o1', zuo.ports['o1'])
            you.connect('o1', s0.ports['o2'])
            s2.connect('o2', you.ports['o2'])
            length = 2 * 97.66 + 3 * mid_length + length
            c.add_port('o1', center=s2.ports['o1'].center,width=s2.ports['o1'].width,orientation=180,layer=gf.get_layer('PORT'))
            c.add_port('o2', center=s1.ports['o2'].center,width=s1.ports['o2'].width,orientation=0,layer=gf.get_layer('PORT'))
        else:
            youori=gf.Component()
            youref = youori<<gf.path.extrude(r, layer=gf.get_layer('WG'), width=width)
            youori.add_port('o3',center=youref.ports['o1'].center,width=youref.ports['o1'].width,orientation=180,layer=gf.get_layer('PORT'))
            youori.add_port('o4', center=youref.ports['o2'].center, width=youref.ports['o2'].width, orientation=180,
                            layer=gf.get_layer('PORT'))
            you=c<<youori
            zuoori=gf.Component()
            zuoref = zuoori<<gf.path.extrude(l, layer=gf.get_layer('WG'), width=width)
            zuoori.add_port('o3',center=zuoref.ports['o1'].center,width=zuoref.ports['o1'].width,orientation=0,layer=gf.get_layer('PORT'))
            zuoori.add_port('o4',center=zuoref.ports['o2'].center,width=zuoref.ports['o2'].width,orientation=0,layer=gf.get_layer('PORT'))
            zuo=c<<zuoori

            zuo.connect('o4', c.ports['o{}'.format((2 * i - 3))])

            you.connect('o3', c.ports['o{}'.format((2 * i - 2))])
            s1 = c << gf.components.straight(width=width, length=mid_length)
            s2 = c << gf.components.straight(width=width, length=mid_length)
            # s1.connect('o1', zuo.ports['o3'])
            # s2.connect('o2', you.ports['o4'])
            s1.center=zuo.ports['o3'].center+[mid_length/2,0]
            s2.center = you.ports['o4'].center + [-mid_length / 2, 0]
            length = length + 2 * mid_length + 8 * L0
            c.add_port('o{}'.format((2 * i - 1)), center=s2.ports['o1'].center, orientation=180,
                       width=s2.ports['o1'].width, layer=gf.get_layer('PORT'))
            c.add_port('o{}'.format((2 * i)), center=s1.ports['o2'].center, orientation=0, width=s1.ports['o2'].width,
                       layer=gf.get_layer('PORT'))

        c.absorb(zuo)
        c.absorb(you)
        c.absorb(s1)
        c.absorb(s2)
    p=gf.Component()
    ref=p<<c
    p.add_port('o1', port=s2.ports['o1'])
    p.add_port('o2', port=s1.ports['o2'])
    print(length)

    return p