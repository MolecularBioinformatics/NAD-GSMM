#!/usr/bin/env python
# coding: utf-8

import svgutils as sg
from pathlib import Path
a4 = ('210mm', '297mm')
bmc_full = ('170', '225')
bmc_half = ('85mm', '225mm')

supplement_path = Path('./images/supplements')

# Fraction of Optimum NAD
atp = sg.compose.SVG(supplement_path / 'fraction_of_optimum_nad_atp.svg')
ldh = sg.compose.SVG(supplement_path / 'fraction_of_optimum_nad_ldh.svg')
c5 = sg.compose.SVG(supplement_path / 'fraction_of_optimum_nad_c5.svg')

x1 = 2
x2 = 80
y1 = 5
y2 = 55
text_size = 5.
scale = .19
offset_x = -4
offset_y = -10

fig = sg.compose.Figure(
    *('160', '105'),
    sg.compose.Panel(
        atp.scale(scale).move(offset_x, offset_y),
        sg.compose.Text('(a)', 0, 0, size=text_size)
    ).move(x1, y1), 
    sg.compose.Panel(
        ldh.scale(scale).move(offset_x, offset_y),
        sg.compose.Text('(b)', 0, 0, size=text_size)
    ).move(x2, y1),
    sg.compose.Panel(
        c5.scale(scale).move(offset_x, offset_y),
        sg.compose.Text('(c)', 0, 0, size=text_size)
    ).move(x1, y2),
)
fig.save(supplement_path / 'figure_s5.svg')

# fraction of optimum GIMME

atp = sg.compose.SVG(supplement_path / 'atp_flux.svg')
n_rxns = sg.compose.SVG(supplement_path / 'number_reactions.svg')

x1 = 5
x2 = 95
y1 = 2
y2 = 57
offset_x = -4
offset_y = 5

scale1 = .19
scale2 = .2
text_size = 5.

fig = sg.compose.Figure(
    *('185', '70'),
    sg.compose.Panel(
        atp.scale(scale1).move(0,0),
        sg.compose.Text('(a)', offset_x, offset_y, size=text_size)
    ).move(x1, y1),
    sg.compose.Panel(
        n_rxns.scale(scale1).move(0, 0),
        sg.compose.Text('(b)', offset_x, offset_y, size=text_size)
    ).move(x2, y1),
)
fig.save(supplement_path / "figure_s3.svg")


# KMs
b = sg.compose.SVG(supplement_path / 'kms_brenda.svg')
s = sg.compose.SVG(supplement_path / 'kms_sabiork.svg')

x1 = 0
x2 = 90
y1 = 0
text_size = 5.
scale = .2
offset_y = 5
offset_x = 2

fig = sg.compose.Figure(
    *('180', '65'),
    sg.compose.Panel(
        b.scale(scale).move(0,0),
        sg.compose.Text('(a)', offset_x, offset_y, size=text_size)
    ).move(x1, y1), 
    sg.compose.Panel(
        s.scale(scale).move(0,0),
        sg.compose.Text('(b)', offset_x, offset_y, size=text_size)
    ).move(x2, y1),
)
fig.save(supplement_path / "figure_s6.svg")
