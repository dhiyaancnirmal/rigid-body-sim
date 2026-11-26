"""Render the demo scene to demo/demo.gif.

A chain of hinged links with a heavy ball swings into a pyramid of boxes
while circles and triangles rain down. Runs at a fixed 60 Hz step and
prints the final state hash, so two runs can be compared bit for bit.

    PYTHONPATH=python python demo/make_gif.py [--frames-dir DIR]
"""

import argparse
import math
import os
import sys

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import rbsim  # noqa: E402

W, H = 640, 400
SCALE = 26.0  # pixels per metre
ORIGIN = (W * 0.42, H - 30)  # world (0, 0) in pixels
DT = 1.0 / 60.0
SECONDS = 9.0

BG = (250, 249, 246)
INK = (40, 40, 48)
COLORS = {
    "static": (120, 120, 128),
    "box": (226, 138, 76),
    "circle": (84, 140, 214),
    "tri": (102, 178, 120),
    "chain": (90, 90, 104),
    "ball": (196, 72, 84),
}


def build():
    w = rbsim.World()
    kinds = []

    def add(kind, idx):
        kinds.append(kind)
        return idx

    add("static", w.add_box(0.0, -0.5, 30.0, 0.5, static=True, friction=0.6))

    # Pyramid of boxes, 8 rows.
    rows, size = 8, 0.5
    for r in range(rows):
        for c in range(rows - r):
            x = 1.0 + (c + 0.5 * r) * 2.0 * size * 1.02
            y = size + r * 2.0 * size
            add("box", w.add_box(x, y, size, size, friction=0.6))

    # Chain: static pivot, 8 links of 0.8 m, then a heavy ball.
    pivot = (-1.5, 10.0)
    prev = add("static", w.add_circle(pivot[0], pivot[1], 0.12, static=True))
    link_len, n_links = 0.8, 8
    x = pivot[0]
    for _ in range(n_links):
        cx = x - 0.5 * link_len
        link = add("chain", w.add_box(cx, pivot[1], 0.5 * link_len, 0.08, density=2.0, friction=0.3))
        w.add_revolute_joint(prev, link, (x, pivot[1]))
        prev = link
        x -= link_len
    r_ball = 0.7
    ball = add("ball", w.add_circle(x - r_ball, pivot[1], r_ball, density=6.0, friction=0.3))
    w.add_revolute_joint(prev, ball, (x, pivot[1]))

    # Rain: circles and triangles above the pyramid.
    for i in range(24):
        bx = 0.5 + (i % 8) * 1.1 + 0.3 * (i // 8)
        by = 12.0 + (i // 8) * 1.4
        if i % 2 == 0:
            add("circle", w.add_circle(bx, by, 0.35, restitution=0.3))
        else:
            add("tri", w.add_polygon([(bx - 0.4, by - 0.3), (bx + 0.4, by - 0.3), (bx, by + 0.4)]))
    return w, kinds


def to_px(x, y):
    return (ORIGIN[0] + SCALE * x, ORIGIN[1] - SCALE * y)


def draw(w, kinds, t, img):
    d = ImageDraw.Draw(img)
    d.rectangle([0, 0, W, H], fill=BG)
    pos = w.positions()
    rot = w.rotations()
    for i, kind in enumerate(kinds):
        x, y = pos[i]
        c, s = rot[i]
        shape = w.shape(i)
        fill = COLORS[kind]
        if shape[0] == "circle":
            r = shape[1]
            px, py = to_px(x, y)
            rp = r * SCALE
            d.ellipse([px - rp, py - rp, px + rp, py + rp], fill=fill, outline=INK, width=1)
            # A radius line shows the spin.
            d.line([to_px(x, y), to_px(x + c * r, y + s * r)], fill=INK, width=1)
        else:
            pts = [to_px(x + c * vx - s * vy, y + s * vx + c * vy) for vx, vy in shape[1]]
            d.polygon(pts, fill=fill, outline=INK)
    d.text((10, 8), f"t = {t:4.2f} s   {len(kinds)} bodies   60 Hz fixed step", fill=INK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "demo.gif"))
    ap.add_argument("--frames-dir", default=None, help="also save every 30th frame as PNG")
    args = ap.parse_args()

    w, kinds = build()
    frames = []
    steps = int(SECONDS / DT)
    for k in range(steps + 1):
        if k % 2 == 0:  # 30 fps
            img = Image.new("RGB", (W, H))
            draw(w, kinds, k * DT, img)
            frames.append(img.quantize(colors=32, method=Image.Quantize.MEDIANCUT))
            if args.frames_dir and k % 60 == 0:
                os.makedirs(args.frames_dir, exist_ok=True)
                img.save(os.path.join(args.frames_dir, f"frame_{k:04d}.png"))
        if k < steps:
            w.step(DT)
    frames[0].save(args.out, save_all=True, append_images=frames[1:], duration=33, loop=0,
                   optimize=True)
    print(f"wrote {args.out}: {len(frames)} frames, final state hash 0x{w.state_hash():016x}")


if __name__ == "__main__":
    main()
