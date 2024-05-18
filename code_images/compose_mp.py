#!/usr/bin/env python
# coding: utf-8
from pathlib import Path
import svgutils as sg


base_path = Path("./images")
img_path = base_path / "parp"
ill_path = base_path
supp_path = base_path / "supplements"
rxns_path = img_path / "rxns_per_path"

for p in (img_path, ill_path, supp_path, rxns_path):
    try:
        assert p.exists()
    except AssertionError as exc:
        print(f"Missing: {p}")
        raise exc

rxns_images = list(rxns_path.glob("*.svg"))
rxn_images = {img.stem: img for img in rxns_images}


# pathway reactions show in manuscript
img_size = 0.075
text_size = 5.0
scale_paths = 0.06


x1 = 0
x2 = 55
y1 = 0
y2 = 65
y3 = 150

fig = sg.compose.Figure(
    *("110mm", "135mm"),
    # b) path selection
    sg.compose.Panel(
        sg.compose.SVG(img_path / "pathways_top_select_rel.svg")
        .scale(scale_paths)
        .move(0, 5),
        sg.compose.Text("(a)", 1, 5, size=text_size),
    ).move(x1, y1),
    # a) Glycolysis
    sg.compose.Panel(
        sg.compose.SVG(rxn_images["glycolysis_rel"]).scale(img_size).move(0, 5),
        sg.compose.Text("(b)", 1, 5, size=text_size),
    ).move(x2, y1),
    # b) TCA
    sg.compose.Panel(
        sg.compose.SVG(rxn_images["tca_cycle_rel"]).scale(img_size).move(0, 5),
        sg.compose.Text("(c)", 1, 5, size=text_size),
    ).move(x1, y2),
    # c) Mal-asp
    sg.compose.Panel(
        sg.compose.SVG(rxn_images["mal_asp_shuttle_rel"]).scale(img_size).move(0, 5),
        sg.compose.Text("(d)", 1, 5, size=text_size),
    ).move(x2, y2),
)
fig.save(ill_path / "figure_2.svg")

# absolute numbers for manuscript pathways to go into supplements
img_size = 0.075
text_size = 5.0

x1 = 0
x2 = 55
y1 = 0
y2 = 60

fig = sg.compose.Figure(
    *("120mm", "105mm"),
    # a) Glycolysis
    sg.compose.Panel(
        sg.compose.SVG(rxn_images["glycolysis_abs"]).scale(img_size).move(0, 5),
        sg.compose.Text("(a)", 1, 5, size=text_size),
    ).move(x1, y1),
    # b) TCA
    sg.compose.Panel(
        sg.compose.SVG(rxn_images["tca_cycle_abs"]).scale(img_size).move(0, 5),
        sg.compose.Text("(b)", 1, 5, size=text_size),
    ).move(x2, y1),
    # c) Mal-asp
    sg.compose.Panel(
        sg.compose.SVG(rxn_images["mal_asp_shuttle_abs"]).scale(img_size).move(0, 5),
        sg.compose.Text("(c)", 1, 5, size=text_size),
    ).move(x1, y2),
)
fig.save(supp_path / "figure_s2.svg")
