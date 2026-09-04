"""Integrate ply stresses from a ccx *EL PRINT S block into a section moment.

Case-local and hard-coded for this model (4 x 0.25 mm plies, 20 mm width, the
ori_p00/ori_p90 orientation names). It is here because ccx will not report a
reaction moment any other way for a composite shell, and because the two things
it gets right -- the integration-point layout and the non-co-rotating stress
frame -- were expensive to work out.

Not promoted into src/compfea yet: the fin's load case is rotation-driven too,
so this belongs in the package as a tested total_moment(), with the energy
cross-check as its regression test. Verified against cantilever_89deg.inp:
mid-span M = -3098.5 N.mm, against 2U/theta = 3098.3 from the energy output.

Usage: python moment_probe.py job.dat
"""
import re, sys, math

A = 1.0/math.sqrt(3.0)

def parse_last_stress_block(path, elset="BLADE"):
    """Return (time, {elem: [(ip, sxx,syy,szz,sxy,sxz,syz, oriname), ...]})"""
    blocks = []
    cur = None
    hdr = re.compile(r"stresses \(elem, integ\.pnt\.,.*\) for set (\S+) and time\s+(\S+)")
    with open(path) as f:
        for line in f:
            m = hdr.search(line)
            if m:
                if m.group(1) == elset:
                    cur = (float(m.group(2)), {})
                    blocks.append(cur)
                else:
                    cur = None
                continue
            if cur is None:
                continue
            t = line.split()
            if len(t) == 9 and t[0].isdigit() and t[1].isdigit():
                e = int(t[0]); ip = int(t[1])
                vals = [float(x) for x in t[2:8]]
                cur[1].setdefault(e, []).append((ip, vals, t[8]))
            elif line.strip() and not line.startswith(" " * 4):
                cur = None
    return blocks[-1]

def ply_z_centres(thicknesses):
    tot = sum(thicknesses)
    z = -tot / 2.0
    out = []
    for t in thicknesses:
        out.append((z + t / 2.0, t))
        z += t
    return out

def section_moment(ips, thicknesses, width, axial_of_ori):
    """ips: list of (ip, [sxx..syz], oriname) for one shell element, 8*nply entries.
    Returns (M_eta_minus, M_eta_plus, M_mean) about the width axis."""
    nply = len(thicknesses)
    d = {ip: (v, o) for ip, v, o in ips}
    assert len(d) == 8 * nply, "expected %d ips, got %d" % (8 * nply, len(d))
    zc = ply_z_centres(thicknesses)
    # within a ply block of 8: local hex ordering
    #   xi  (width) : pts 1,4,5,8 = -a ; 2,3,6,7 = +a
    #   eta (length): pts 1,2,5,6 = -a ; 3,4,7,8 = +a
    #   zeta(thick) : pts 1..4     = -a ; 5..8    = +a
    ETA = {1: 0, 2: 0, 3: 1, 4: 1, 5: 0, 6: 0, 7: 1, 8: 1}
    ZET = {1: -A, 2: -A, 3: -A, 4: -A, 5: A, 6: A, 7: A, 8: A}
    M = [0.0, 0.0]
    for p in range(nply):
        zcen, t = zc[p]
        for k in range(1, 9):
            ip = 8 * p + k
            vals, ori = d[ip]
            sig = axial_of_ori(ori, vals)
            z = zcen + (t / 2.0) * ZET[ip - 8 * p]
            # weights 1*1, jacobian (width/2)*(t/2)
            M[ETA[k]] += sig * z * (width / 2.0) * (t / 2.0)
    return M[0], M[1], 0.5 * (M[0] + M[1])

def axial_cfrp(ori, vals):
    """Axial (beam-direction) Cauchy stress.

    ccx prints stresses in the *ORIENTATION frame, which does NOT co-rotate with
    the material under NLGEOM.  The beam bends about global x, so for both plies
    the beam axis stays in the plane spanned by (beam-axis-at-t0, shell-normal),
    and the trace of that 2x2 block is invariant under the rotation:
        s_axial + s_shellnormal, with s_shellnormal ~ 0 in a shell.
    ori_p00: local1=+y(beam) local2=-x local3=+z(normal)  ->  sxx + szz
    ori_p90: local1=+x local2=+y(beam) local3=+z(normal)  ->  syy + szz
    """
    u = ori.upper()
    if "P00" in u:
        return vals[0] + vals[2]     # sxx + szz
    if "P90" in u:
        return vals[1] + vals[2]     # syy + szz
    raise ValueError(ori)

if __name__ == "__main__":
    path = sys.argv[1]
    th = [0.25] * 4
    time, els = parse_last_stress_block(path)
    print("last stress block at time %g" % time)
    for e in sorted(els):
        m0, m1, mm = section_moment(els[e], th, 20.0, axial_cfrp)
        print("  elem %3d   M(eta-)=%10.3f  M(eta+)=%10.3f  M(mean)=%10.3f" % (e, m0, m1, mm))
