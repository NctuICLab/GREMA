#!/usr/bin/env python3
"""Make the rows of an external input for a GREMA profile (GREMA option -input) from a small table, and append them to the profile.

usage: make_input.py DATA.txt TABLE.tsv NAME [--times t0,t1,...] [--shape step|decay] [--halftime H] [--reverse] [--out NEW.txt]
  DATA.txt  the GREMA profile (repeat_number=, timepoint_number=, rows 'rep<k> <gene> v1 ... vT'); it gives the number of series and time points
  TABLE.tsv one line per series, tab separated, with a header line:  series  start  level  [baseline]
              series    1, 2, ... (every series of DATA.txt must have a line)
              start     when the stimulus begins: the index of the time point (0 is the first) or, with --times, a time in the unit of --times
              level     the raw strength after the start (a dose, a concentration, 1 for 'on')
              baseline  the raw strength before the start, default 0
            The raw values are scaled with one linear map through zero to the scale of the genes: u = 0.01 + 0.99 * raw / (the largest raw value of the table),
            so that the same raw value gets the same u in every series and 0 becomes 0.01 (GREMA cannot use 0).
  --times   the real sampling times of the time points (the list of the GREMA option -time); then 'start' is a time and the stimulus is on from the first time
            point at or after it; the step of the model between two points stays 1 (as in GREMA)
  --shape   step (default): the level from the start on; decay: the level at the start, then falling to the baseline with the half-time --halftime
            (in the unit of --times, or in time points) like  baseline + (level - baseline) * 2^(-(t - start) / halftime)
  --reverse the values of every series in the reverse order of time (a control for the shape: a stimulus that grows instead of decaying)
  --out     the profile with the rows of NAME added (default: NAME rows only, to the standard output)
The data check (check_data.py --input NAME) tells whether the input can be identified; at least two series with a different level or start time are needed.
Example (two series, a stimulus at time point 4 with level 1 and 0.5):  see evaluation/analysis/examples/irma_input/README.md
"""
import sys, math


def main():
    a = sys.argv[1:]
    opt = {"--times": None, "--shape": "step", "--halftime": None, "--out": None, "--reverse": None}
    pos = []
    i = 0
    while i < len(a):
        if a[i] == "--reverse":
            opt["--reverse"] = True; i += 1
        elif a[i] in opt:
            opt[a[i]] = a[i + 1]; i += 2
        else:
            pos.append(a[i]); i += 1
    if len(pos) != 3:
        print(__doc__); sys.exit(1)
    data, table, name = pos
    nrep = npt = None
    for l in open(data):
        if l.startswith("repeat_number="):
            nrep = int(l.split("=")[1])
        elif l.startswith("timepoint_number="):
            npt = int(l.split("=")[1])
    if nrep is None or npt is None:
        sys.exit("%s has no repeat_number= / timepoint_number= lines" % data)
    times = [float(x) for x in opt["--times"].split(",")] if opt["--times"] else list(range(npt))
    if len(times) != npt:
        sys.exit("--times has %d values, the profile has %d time points" % (len(times), npt))
    if opt["--shape"] not in ("step", "decay"):
        sys.exit("--shape must be step or decay")
    if opt["--shape"] == "decay" and not opt["--halftime"]:
        sys.exit("--shape decay needs --halftime")
    rows = {}
    for l in open(table):
        f = l.rstrip("\n").rstrip("\r").split("\t")
        if not f or f[0] in ("", "series") or f[0].startswith("#"):
            continue
        if len(f) < 3:
            sys.exit("%s: expected series, start, level [, baseline]: %r" % (table, l))
        rows[int(f[0])] = (float(f[1]), float(f[2]), float(f[3]) if len(f) > 3 and f[3] != "" else 0.0)
    missing = [r for r in range(1, nrep + 1) if r not in rows]
    if missing:
        sys.exit("%s has no line for the series %s (the profile has %d series)" % (table, missing, nrep))
    extra = [r for r in rows if r > nrep or r < 1]
    if extra:
        sys.exit("%s has lines for the series %s that are not in the profile" % (table, extra))
    raw_max = max(max(abs(l), abs(b)) for _, l, b in rows.values())
    if raw_max <= 0:
        sys.exit("all levels are 0")
    half = float(opt["--halftime"]) if opt["--halftime"] else None
    out = []
    for r in range(1, nrep + 1):
        start, level, base = rows[r]
        v = []
        for t in times:
            if t < start:
                raw = base
            elif opt["--shape"] == "step":
                raw = level
            else:
                raw = base + (level - base) * 2.0 ** (-(t - start) / half)
            v.append(0.01 + 0.99 * raw / raw_max)
        if opt["--reverse"]:
            v = v[::-1]  # a control: the same values in the reverse order of time
        out.append("rep%d\t%s\t%s" % (r, name, "\t".join("%.6f" % x for x in v)))
    if opt["--out"]:
        body = open(data).read().rstrip("\n").split("\n")
        if any(len(l.split("\t")) > 2 and l.split("\t")[1] == name for l in body):
            sys.exit("%s already has rows named %s" % (data, name))
        with open(opt["--out"], "w") as f:
            f.write("\n".join(body + out) + "\n")
        print("wrote %s: %d series, %d time points, input %s (raw max %g, shape %s)" % (opt["--out"], nrep, npt, name, raw_max, opt["--shape"]))
    else:
        print("\n".join(out))


main()
