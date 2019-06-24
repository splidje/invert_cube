import sys
import re
import numpy
from collections import defaultdict

in_path = sys.argv[1]
out_path = sys.argv[2]

in_file = open(in_path)

size = None
while size is None:
    line = in_file.readline()
    match = re.match(r"\s*LUT_3D_SIZE\s+(\d+)", line)
    if match:
        size = int(match.group(1))

inverse_map = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

for b in range(size):
    for g in range(size):
        for r in range(size):
            R = None
            while R is None:
                line = in_file.readline()
                match = re.match(r"\s*([+-.\d]+)\s+([+-.\d]+)\s+([+-.\d]+)", line)
                if match:
                    R, G, B = tuple(int(float(x) * (size-1)) for x in match.groups())
            inverse_map[R][G][B].append((r, g, b))
in_file.close()

# average
for r, gbs in inverse_map.iteritems():
    for g, bs in gbs.iteritems():
        for b, vals in bs.iteritems():
            val = tuple(sum(x)/len(x) for x in zip(*vals))
            inverse_map[r][g][b] = val

# set up ocio
import PyOpenColorIO as ocio
conf = ocio.Config()
conf.addColorSpace(ocio.ColorSpace('ref'))
ft = ocio.FileTransform(in_path)
ft.setInterpolation('tetrahedral')
cs = ocio.ColorSpace('test')
cs.setTransform(ft, ocio.Constants.COLORSPACE_DIR_FROM_REFERENCE)
conf.addColorSpace(cs)
proc = conf.getProcessor('ref', 'test')

out_file = open(out_path, "w")

out_file.write("LUT_3D_SIZE {}\n".format(size))

tris = {}

for b in range(size):
    for g in range(size):
        for r in range(size):
            R = G = B = 0
            val = inverse_map[r][g][b]
            if len(val) == 0:
               
                # INTERPOLATION!!
                
                # find a tetrahedron we're in
                minmaxs = ((x, x) for x in (r, g, b))
                points = []
                tetras = set()
                tetra = None
                while not tetra:
                    # expand
                    minmaxs = tuple(
                        (
                            (mn - 1 if mn > 0 else mn),
                            (mx + 1 if mx < (size-1) else mx),
                        ) for mn, mx in minmaxs
                    )
                    # add points
                    for rr in minmaxs[0]:
                        for gg in minmaxs[1]:
                            for bb in minmaxs[2]:
                                if len(inverse_map[rr][gg][bb]) > 0:
                                    points.append((rr, gg, bb))
                    points = list(set(points))
                    # find all tetras
                    for i, p1 in enumerate(points):
                        for j, p2 in enumerate(points[i+1:]):
                            for k, p3 in enumerate(points[i+j+2:]):
                                for l, p4 in enumerate(points[i+j+k+3:]):
                                    key = tuple(sorted((p1, p2, p3, p4)))
                                    if key not in tetras:
                                        # middle
                                        cntr = [sum(p)/float(len(p)) for p in zip(*key)]
                                        # establish faces
                                        faces = []
                                        for ii, t1 in enumerate(key):
                                            for jj, t2 in enumerate(key[ii+1:]):
                                                for kk, t3 in enumerate(key[ii+jj+2:]):
                                                    faces.append(tuple(sorted((t1, t2, t3))))
                                        # right side of each face?
                                        inside = True
                                        for face in faces:
                                            # plane eq
                                            if face not in tris:
                                                # get normal
                                                v1 = numpy.subtract(face[1], face[0])
                                                v2 = numpy.subtract(face[2], face[0])
                                                tris[face] = numpy.cross(v1, v2)
                                            norm = tris[face]
                                            # which side centre?
                                            cdist = numpy.dot(numpy.subtract(cntr, face[0]), norm)
                                            if cdist == 0:
                                                inside = False
                                                break
                                            else:
                                                # which side us?
                                                dist = numpy.dot(numpy.subtract((r, g, b), face[0]), norm)
                                                if cdist < 0 and dist > 0 or cdist > 0 and dist < 0:
                                                    inside = False
                                                    break
                                        if inside:
                                            tetra = key
                                            break
                                        tetras.add(key)
                                    if tetra:
                                        break
                                if tetra:
                                    break
                            if tetra:
                                break
                        if tetra:
                            break
                
                # found tetra!

                rgb = r,g,b

                val = None

                # generate pixels in this range
                ranges = [(min(x), max(x)-min(x)) for x in zip(*tetra)]
                pixels = []
                for k in range(size):
                    for j in range(size):
                        for i in range(size):
                            rr = (i * ranges[0][1]/float(size-1) + ranges[0][0]) / (size-1)
                            gg = (j * ranges[1][1]/float(size-1) + ranges[1][0]) / (size-1)
                            bb = (k * ranges[2][1]/float(size-1) + ranges[2][0]) / (size-1)
                            pixels.extend((rr, gg, bb))
                inv_pixels = proc.applyRGB(pixels)
                for i in range(0, len(inv_pixels), 3):
                    rrggbb = tuple(int(round(x*(size-1))) for x in inv_pixels[i:i+3])
                    if (rrggbb == rgb):
                        print "found!"
                        val = pixels[i:i+3]
                        break

                if val is None:
                    print "not found :-("
                    val = (0,0,0)
                    #sys.exit(-1)
                
                """
                # get mapped tetra
                inv_tetra = tuple(inverse_map[rr][gg][bb] for rr, gg, bb in tetra)
                

                face1 = tetra[:3]
                norm1 = tris[face1]
                norm1 = numpy.true_divide(norm1, numpy.linalg.norm(norm1))
                face1d = numpy.dot(face1[0], norm1)
                p4 = tetra[3]
                # draw a line from point 4 (not on face 1) to our point
                v3 = numpy.subtract(rgb, p4)
                # get where it hits opposite plane (face 1) - get distance from point 4 as proportion of distance to new point
                mult = (face1d - numpy.dot(p4, norm1)) / numpy.dot(v3, norm1)
                if mult < 0.99999:
                    print "fucked1!", mult
                    sys.exit(-1)
                move3 = 1 - 1 / mult # proportional distance from face
                rgb2 = numpy.add(numpy.multiply(v3, mult), p4)
                # draw a line from point 3 on face 1 to this new point
                p3 = face1[2]
                if tuple(map(round, rgb2)) == tuple(p3):
                    rgb1 = p3
                    move2 = 1
                else:
                    v2 = numpy.subtract(rgb2, p3)
                    # get where it hits plane of face made with point 4 and points 1 and 2 on face 1 - get distance from point 3 as proportion of distance to this even newer point
                    face2 = face1[:2] + (p4,)
                    norm2 = tris[face2]
                    norm2 = numpy.true_divide(norm2, numpy.linalg.norm(norm2))
                    face2d = numpy.dot(face2[0], norm2)
                    mult = (face2d - numpy.dot(p3, norm2)) / numpy.dot(v2, norm2)
                    if mult < 0.9999:
                        #print "fucked2!", v2, norm2, mult
                        rgb1 = p3
                        move2 = 1
                        #sys.exit(-1)
                    elif numpy.isinf(mult):
                        print "inf", face2, norm2, face2d, p3, rgb2, v2, tuple(map(round, rgb2)), tuple(p3), tuple(map(round, rgb2)) == tuple(p3)
                        sys.exit(-1)
                    else:
                        move2 = 1 - 1 / mult # proporitional distance from edge
                        rgb1 = numpy.add(numpy.multiply(v2, mult), p3)
                # get distance from this even newer point to point 1 on face 1 - as a proportion of distance to point 2 from point 1
                move1 = numpy.linalg.norm(numpy.subtract(rgb1, face1[0])) / numpy.linalg.norm(numpy.subtract(face1[1], face1[0]))

                # follow these three moves on the inverse tetra
                # move from p1 on f1 towards p2
                face1 = inv_tetra[:3]
                # get normal
                v1 = numpy.subtract(face1[1], face1[0])
                v2 = numpy.subtract(face1[2], face1[0])
                norm1 = numpy.cross(v1, v2)
                norm1d = numpy.linalg.norm(norm1)
                if not numpy.isfinite(norm1d):
                    print norm1, norm1d, v1, v2
                    sys.exit(-1)
                norm1 = numpy.true_divide(norm1, norm1d)
                face1d = numpy.dot(face1[0], norm1)
                p1, p2 = face1[:2]
                v1 = numpy.subtract(p2, p1)
                rgb1 = numpy.add(p1, numpy.multiply(v1, move1))
                # move from rgb1 to p3 on f1
                p3 = face1[2]
                v2 = numpy.subtract(p3, rgb1)
                rgb2 = numpy.add(rgb1, numpy.multiply(v2, move2))
                # move from rgb2 to p4 on tetra
                p4 = inv_tetra[3]
                v3 = numpy.subtract(p4, rgb2)
                rgb3 = numpy.add(rgb2, numpy.multiply(v3, move3))

                val = rgb3
                """
                
            R, G, B = tuple((x / float(size-1)) for x in val)
            out_file.write("{} {} {}\n".format(R, G, B))
out_file.close()

