#!/usr/bin/env python
import os
import sys
import re
from collections import defaultdict
import PyOpenColorIO as ocio
import numpy
import time
import itertools
import subprocess

in_cube_path = sys.argv[1]
out_cube_path = sys.argv[2]

conf = ocio.Config()
conf.addColorSpace(ocio.ColorSpace('ref'))
ft = ocio.FileTransform(in_cube_path)
ft.setInterpolation('tetrahedral')
cs = ocio.ColorSpace('test')
cs.setTransform(ft, ocio.Constants.COLORSPACE_DIR_FROM_REFERENCE)
conf.addColorSpace(cs)
proc = conf.getProcessor('ref', 'test')

# get cube size
in_cube_f = open(in_cube_path)
size = None
while size is None:
    line = in_cube_f.readline()
    match = re.match(r"^\s*LUT_3D_SIZE\s+(\d+)", line)
    if match:
        size = int(match.group(1))
in_cube_f.close()

# == FUNCTIONS ==
def transform_range(ranges, precision):
    pixels = []
    sizes = [(mx - mn) for mn, mx in ranges]
    for b in range(precision):
        for g in range(precision):
            for r in range(precision):
                pixels.extend((
                    (x * sizes[i]) / float(precision-1) + ranges[i][0]
                    for i, x in enumerate((r, g, b))
                ))
    trans_pix = proc.applyRGB(pixels)
    ret = []
    for i in range(0, len(trans_pix), 3):
        ret.append((tuple(pixels[i:i+3]), tuple(trans_pix[i:i+3])))
    return ret

# == MAIN STUFF ==

st = time.time()

# do first transform
val_table = transform_range([(0, 1)] * 3, size)

# get output min and maxes
ranges = [(min(x), max(x)) for x in zip(*(d for s, d in val_table))]
#ranges = [(0, 1)] * 3

# TODO: flip this float table, and collect lists of multiple destinations
inv_val_map = defaultdict(list)
for s, d in val_table:
    inv_val_map[d].append(s)
# TODO: connect up smallest unique tetras source side

pnts = inv_val_map.keys()

# write out points
path_prefix = os.path.splitext(out_cube_path)[0]
nodes_path = "{}.node".format(path_prefix)
node_f = open(nodes_path, "w")
node_f.write("{} 3 0 0\n".format(len(pnts)))
for i, d in enumerate(pnts):
    node_f.write("{} {}\n".format(i + 1, " ".join(map(str, d))))
node_f.close()
# run tetgen
os.system(subprocess.list2cmdline(('tetgen', '-M', nodes_path)))
# read in tetras
tetras_path = "{}.1.ele".format(path_prefix)
tetras_f = open(tetras_path)
num_tetras = int(re.search(r"^\s*(\d+)", tetras_f.readline()).group(1))
tetras = tuple(tuple(map(int, re.findall(r"\d+", tetras_f.readline())[1:5])) for i in range(num_tetras))
tetras_f.close()

print time.time() - st

# TODO: scale up and round dests to find gaps
inv_map_scaled = defaultdict(list)
for d in pnts:
    inv_map_scaled[tuple(
        int(round((size-1) * (x - r[0]) / (r[1] - r[0])))
        for x, r in zip(d, ranges)
    )].extend(inv_val_map[d])
# aggregate
for RGB, ss in inv_map_scaled.items():
    inv_map_scaled[RGB] = [numpy.mean(x) for x in zip(*ss)]
        #sorted(ss)[0]
        #ss[0] if len(ss) == 1 else (1, 0.5, 0.5)

# IDEA: DO EVERYTHING WITH TETRAHEDRONS!! AHAHAHAHAHAHA! USE THE SAME PRINCIPLE FOR EVERY SINGLE M.F. POINT WE WANT TO MAP!

gaps = []
for R in range(size):
    for G in range(size):
        for B in range(size):
            RGB = R, G, B
            #if not inv_map_scaled[RGB]
            gaps.append(tuple(r[0] + (r[1] - r[0]) * x / float(size-1) for x, r in zip(RGB, ranges)))

print time.time() - st

def in_tetra(tetra_dets, p1):
    for face0, nrm, pdpos in tetra_dets:
        p1d = numpy.dot(nrm, numpy.subtract(p1, face0))
        if p1d != 0 and (p1d > 0) != pdpos:
            return False
    return True

tetra_cache = {}

def get_tetra_dets(idx):
    if idx not in tetra_cache:
        tetra = tuple(pnts[i-1] for i in tetras[idx])
        tetra_cache[idx] = tetra, []
        for i in range(4):
            face = [tetra[j] for j in range(4) if j != i]
            nrm = numpy.cross(numpy.subtract(face[1], face[0]), numpy.subtract(face[2], face[0]))
            pdpos = numpy.dot(nrm, numpy.subtract(tetra[i], face[0])) > 0
            tetra_cache[idx][1].append((face[0], nrm, pdpos))
    return tetra_cache[idx]

def calc_sphere(tetra):
    coeffs = []
    Tdet = numpy.linalg.det([tetra[i] + (1,) for i in range(4)])
    ts = []
    t = [-numpy.dot(tetra[i], tetra[i]) for i in range(4)]
    for i in range(4):
        M = [
            [t[j] if k == i else (tetra[j][k] if k < 3 else 1) for k in range(4)]
            for j in range(4)
        ]
        coeffs.append(numpy.linalg.det(M) / Tdet)
    centre = [-x/2 for x in coeffs[:3]]
    radius = numpy.sqrt(numpy.dot(coeffs[:3], coeffs[:3]) - 4*coeffs[3]) / 2
    return centre, radius

tetra_zones = defaultdict(list)
for idx in range(len(tetras)):
    tetra = tuple(pnts[i-1] for i in tetras[idx])
    # range of grid cells points are in
    gr = [(min(x), max(x) + 1) for x in zip(*((int((size-1) * (x - r[0]) / (r[1] - r[0])) for x, r in zip(p, ranges)) for p in tetra))]
    for r in range(*gr[0]):
        for g in range(*gr[1]):
            for b in range(*gr[2]):
                tetra_zones[(r, g, b)].append(idx)

# TODO: find the tetras each (unscaled float) gap is in:
for gap in gaps:
    inside_tetra = None
    # which grid bit?
    gcell = tuple(int((size-1) * (x - r[0]) / (r[1] - r[0])) for x, r in zip(gap, ranges))
    # which tetra?
    cnds = set()
    for tetra_idx in tetra_zones[gcell]:
        tetra, tetra_dets = get_tetra_dets(tetra_idx)
        # add points to candidates for later
        cnds.update(tetra)
        is_in = in_tetra(tetra_dets, gap)
        if is_in:
            inside_tetra = tetra, tetra_dets
            break
    # found tetra?
    if inside_tetra:
        tetra, tetra_dets = inside_tetra

        # git distance from point 4
        v3 = numpy.subtract(gap, tetra[3])
        face4_dets = tetra_dets[3]
        nrm = face4_dets[1]
        f4d = numpy.dot(nrm, numpy.subtract(face4_dets[0], tetra[3]))
        v3d = numpy.dot(nrm, v3)
        move3 = (v3d / f4d)
        if move3 != 0:
            p2 = numpy.add(tetra[3], numpy.true_divide(v3, move3))
        else:
            p2 = tetra[2]
        move3 = 1 - move3 # reverse move is rest of the way
        # distance from point3
        v2 = numpy.subtract(p2, tetra[2])
        face3_dets = tetra_dets[2]
        nrm = face3_dets[1]
        f3d = numpy.dot(nrm, numpy.subtract(face3_dets[0], tetra[2]))
        v2d = numpy.dot(nrm, v2)
        move2 = (v2d / f3d)
        if move2 != 0:
            p1 = numpy.add(tetra[2], numpy.true_divide(v2, move2))
        else:
            p1 = tetra[0]
        move2 = 1 - move2
        # distance from point1
        move1 = numpy.linalg.norm(numpy.subtract(p1, tetra[0])) / numpy.linalg.norm(numpy.subtract(tetra[1], tetra[0]))

        # map tetra
        # TODO: should really pick smallest
        min_rad = None
        for p1 in inv_val_map[tetra[0]]:
            for p2 in inv_val_map[tetra[1]]:
                for p3 in inv_val_map[tetra[2]]:
                    for p4 in inv_val_map[tetra[3]]:
                        tet = (p1, p2, p3, p4)
                        cntr, radius = calc_sphere(tet)
                        if min_rad is None or radius < min_rad:
                            min_rad = radius
                            tetra_src = tet

        # do the moves
        p1 = numpy.add(tetra_src[0], numpy.multiply(numpy.subtract(tetra_src[1], tetra_src[0]), move1))
        p2 = numpy.add(p1, numpy.multiply(numpy.subtract(tetra_src[2], p1), move2))
        p3 = numpy.add(p2, numpy.multiply(numpy.subtract(tetra_src[3], p2), move3))
        
        rgb = p3

        #rgb = [numpy.mean(x) for x in zip(*tetra_src_pnts)]

        # TODO:     work out weights for each vertex (sum of vertices / 2 is midpoint, sum of 0.5 * each vertex is same, weight per vertex describes everywhere within tetra?)
        # TODO:     every combination of possible mapping vertices forms a tetra. find smallest
        # TODO:     apply the weights to that smallest tetra to get coords
    # didn't find tetra?
    else:
        # closest point
        srted = sorted((numpy.linalg.norm(numpy.subtract(cnd, gap)), cnd) for cnd in cnds)
        if not srted:
            srted = sorted((numpy.linalg.norm(numpy.subtract(cnd, gap)), cnd) for cnd in pnts)
        RGB = srted[0][1]
        # TODO: if nothing in this cell, search around the neighbours            
        # map
        rgb = tuple(numpy.mean(x) for x in zip(*inv_val_map[RGB]))

    inv_map_scaled[
        tuple(
            int(round((size-1) * (x - r[0]) / (r[1] - r[0])))
            for x, r in zip(gap, ranges)
        )
    ] = rgb

print time.time() - st

# write out
out_cube_f = open(out_cube_path, "w")
out_cube_f.write("LUT_3D_SIZE {}\n".format(size))
out_cube_f.write("DOMAIN_MIN {}\n".format(" ".join(str(r[0]) for r in ranges)))
out_cube_f.write("DOMAIN_MAX {}\n".format(" ".join(str(r[1]) for r in ranges)))
for B in range(size):
    for G in range(size):
        for R in range(size):
            RGB = (R, G, B)
            rgb = inv_map_scaled[RGB]
            #if not rgb:
            #    rgb = (0, 0, 0)
            out_cube_f.write("{}\n".format(" ".join(map(str, rgb))))
out_cube_f.close()

