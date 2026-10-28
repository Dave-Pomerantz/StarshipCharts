"""Builds one chart set: NOAA's charts and depth grids for the region, with only the changed files new.

A set is one GitHub release. Its set.tsv lists every file of the set; a file unchanged since an earlier set is
not uploaded again, and its line names the release that holds it. A chart is NOAA's zip as NOAA publishes it
(the base file, all its updates, NOAA's user agreement); it changes when NOAA's edition or update number
changes. A chart NOAA cancels, or that leaves the region, is not in the next set. A depth grid is one NOAA
S-102 file converted by grid_file.py; it changes when NOAA issues the file again under a new name.

  python build_set.py --work DIR --out DIR [--previous URL_OR_PATH] [--tag TAG]

Writes OUT/set.tsv, the new files under OUT/assets, and OUT/upload-<release>.txt per release to create.
Exit code 0 when the set differs from the previous one, 3 when nothing changed, 1 on a failure.
"""
import argparse
import datetime as dt
import hashlib
import os
import re
import sys
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

import grid_file

CHARTS_CATALOGUE = "https://charts.noaa.gov/ENCs/ENCProdCat.xml"
GRIDS_BASE = "https://noaa-s102-pds.s3.amazonaws.com/ed3.0.0/"
GRIDS_CATALOGUE = GRIDS_BASE + "_CATALOG/CATALOG.XML"
RELEASE_ASSETS = 900  # GitHub allows 1,000 files per release; the rest go to <tag>-2, <tag>-3 and on.
COLUMNS = ["kind", "name", "version", "release", "asset", "bytes", "sha256", "west", "south", "east", "north"]
USER_AGENT = "StarshipCharts/1 (+https://github.com/Dave-Pomerantz/StarshipCharts)"


def log(*a):
    print(*a, flush=True)


def fetch(url, path=None, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as r:
                if path is None:
                    return r.read()
                tmp = path + ".partial"
                with open(tmp, "wb") as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        f.write(b)
                os.replace(tmp, path)
                return path
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
            err = e
        except Exception as e:  # network faults: try again
            err = e
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{url}: {err}")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def overlaps(b, region):
    w, s, e, n = b
    rw, rs, re_, rn = region
    return not (e < rw or w > re_ or n < rs or s > rn)


def read_region(path):
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and not line.startswith("west"):
                return tuple(float(x) for x in line.split("\t")[:4])
    raise ValueError("region.tsv holds no box")


def charts_in(xml, region):
    root = ET.fromstring(xml)
    stamp = root.findtext("Header/dt_valid") or ""
    out = {}
    for cell in root.findall("cell"):
        if cell.findtext("status") != "Active":
            continue
        pts = [(float(v.findtext("long")), float(v.findtext("lat"))) for v in cell.iter("vertex")]
        if not pts:
            continue
        b = (min(p[0] for p in pts), min(p[1] for p in pts), max(p[0] for p in pts), max(p[1] for p in pts))
        if not overlaps(b, region):
            continue
        name = cell.findtext("name")
        version = f"{int(cell.findtext('edtn'))}.{int(cell.findtext('updn'))}"
        out[name] = {"kind": "chart", "name": name, "version": version, "url": cell.findtext("zipfile_location"),
                     "asset": f"{name}-{version.replace('.', '-')}.zip", "bounds": b}
    return stamp, out


NS = {"S100XC": "http://www.iho.int/s100/xc/5.2", "gex": "http://standards.iso.org/iso/19115/-3/gex/1.0",
      "gco": "http://standards.iso.org/iso/19115/-3/gco/1.0"}


def grids_in(xml, region):
    root = ET.fromstring(xml)
    out = {}
    for md in root.iter():
        if not md.tag.endswith("S100_DatasetDiscoveryMetadata"):
            continue
        fname = dsid = None
        box = {}
        for e in md.iter():
            tag = e.tag.split("}")[-1]
            if tag == "fileName" and fname is None:
                fname = (e.text or "").strip()
            elif tag == "datasetID" and dsid is None:
                dsid = (e.text or "").strip()
            elif tag in ("westBoundLongitude", "eastBoundLongitude", "southBoundLatitude", "northBoundLatitude") and tag not in box:
                dec = [c for c in e.iter() if c.tag.endswith("Decimal")]
                if dec:
                    box[tag] = float(dec[0].text)
        if not fname or not fname.endswith(".h5") or len(box) < 4:
            continue
        rel = fname.split("file:../", 1)[-1]
        stem = os.path.basename(rel)[:-3]
        # NOAA's name: 102, the producer US00, the tile, two digits of year, three of version, the edition digit.
        m = re.match(r"^(102US00\w+)(\d\d\w{3}\d)$", stem)
        if not m:
            log("skipped, name not understood:", rel)
            continue
        b = (box["westBoundLongitude"], box["southBoundLatitude"], box["eastBoundLongitude"], box["northBoundLatitude"])
        if not overlaps(b, region):
            continue
        sha = dsid.rsplit(":", 1)[-1] if dsid and "sha256" in dsid else None
        name, version = m.group(1), m.group(2)
        old = out.get(name)
        if old and old["version"] >= version:
            continue
        out[name] = {"kind": "grid", "name": name, "version": version, "url": GRIDS_BASE + rel,
                     "source_sha256": sha, "asset": stem + ".grid", "bounds": b}
    return out


def read_set(text):
    rows = {}
    meta = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            parts = line[1:].strip().split("\t")
            if len(parts) >= 2:
                meta[parts[0]] = parts[1]
            continue
        parts = line.split("\t")
        if parts[0] == "kind":
            continue
        r = dict(zip(COLUMNS, parts))
        rows[(r["kind"], r["name"])] = r
    return meta, rows


def previous_set(source):
    if not source:
        return {}, {}
    try:
        if os.path.exists(source):
            with open(source, encoding="utf-8") as f:
                return read_set(f.read())
        return read_set(fetch(source).decode("utf-8"))
    except urllib.error.HTTPError as e:
        if e.code == 404:
            log("no previous set at", source)
            return {}, {}
        raise


def make_new(item, work, assets):
    """Downloads (and for a grid converts) one item into assets; returns its row values."""
    target = os.path.join(assets, item["asset"])
    if item["kind"] == "chart":
        fetch(item["url"], target)
    else:
        src = os.path.join(work, os.path.basename(item["url"]))
        if not (os.path.exists(src) and item.get("source_sha256") and sha256(src) == item["source_sha256"]):
            fetch(item["url"], src)
        if item.get("source_sha256") and sha256(src) != item["source_sha256"]:
            raise RuntimeError(f"{item['url']}: the file does not match NOAA's catalogue checksum")
        grid_file.convert(src, target)
    return os.path.getsize(target), sha256(target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--previous", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--region", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "region.tsv"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0, help="at most this many new files, for trials")
    a = ap.parse_args()
    now = dt.datetime.now(dt.timezone.utc)
    tag = a.tag or now.strftime("set-%Y%m%d-%H%M")
    region = read_region(a.region)
    os.makedirs(a.work, exist_ok=True)
    assets = os.path.join(a.out, "assets")
    os.makedirs(assets, exist_ok=True)

    stamp, charts = charts_in(fetch(CHARTS_CATALOGUE), region)
    grids = grids_in(fetch(GRIDS_CATALOGUE), region)
    log(f"region {region}: {len(charts)} charts, {len(grids)} grids; NOAA charts catalogue {stamp}")
    _, prev = previous_set(a.previous)

    wanted = {("chart", k): v for k, v in charts.items()}
    wanted.update({("grid", k): v for k, v in grids.items()})
    rows, new = [], []
    for key in sorted(wanted):
        item = wanted[key]
        old = prev.get(key)
        if old and old["version"] == item["version"]:
            rows.append(old)
        else:
            new.append(item)
    removed = sorted(k for k in prev if k not in wanted)
    if a.limit:
        new = new[:a.limit]
    log(f"{len(rows)} unchanged, {len(new)} new or changed, {len(removed)} removed")
    for k in removed:
        log("removed:", k[0], k[1], prev[k]["version"])
    if not new and not removed:
        log("nothing changed")
        return 3

    releases = {}
    for i, item in enumerate(new):
        n = i // RELEASE_ASSETS
        item["release"] = tag if n == 0 else f"{tag}-{n + 1}"
        releases.setdefault(item["release"], []).append(item["asset"])

    started = time.time()
    failed = []
    same = []

    def one(item):
        try:
            size, digest = make_new(item, a.work, assets)
            old = prev.get((item["kind"], item["name"]))
            if old and old["sha256"] == digest:
                # NOAA issued the file again with the same depths, uncertainties and surveys: the earlier copy
                # stands, under NOAA's new version, and nothing is uploaded or downloaded again.
                os.remove(os.path.join(assets, item["asset"]))
                same.append(item["asset"])
                return dict(old, version=item["version"])
            b = item["bounds"]
            return {"kind": item["kind"], "name": item["name"], "version": item["version"], "release": item["release"],
                    "asset": item["asset"], "bytes": str(size), "sha256": digest,
                    "west": f"{b[0]:.6f}", "south": f"{b[1]:.6f}", "east": f"{b[2]:.6f}", "north": f"{b[3]:.6f}"}
        except Exception as e:
            log("failed:", item["kind"], item["name"], e)
            failed.append(item)
            return None

    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for i, row in enumerate(pool.map(one, new), 1):
            if row:
                rows.append(row)
            if i % 25 == 0:
                log(f"{i} of {len(new)} in {time.time() - started:.0f} s")
    if failed:
        log(f"{len(failed)} files failed; no set is written")
        return 1

    rows.sort(key=lambda r: (r["kind"], r["name"]))
    with open(os.path.join(a.out, "set.tsv"), "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# set\t{tag}\n")
        f.write(f"# built\t{now.strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
        f.write(f"# noaa charts catalogue\t{stamp}\n")
        f.write(f"# region\t{' '.join(str(x) for x in region)}\n")
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in COLUMNS) + "\n")
    for name in same:
        for names in releases.values():
            if name in names:
                names.remove(name)
    releases = {rel: names for rel, names in releases.items() if names}
    for rel, names in releases.items():
        with open(os.path.join(a.out, f"upload-{rel}.txt"), "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(names) + "\n")
    total = sum(int(r["bytes"]) for r in rows)
    added = sum(int(r["bytes"]) for r in rows if r["release"] in releases)
    log(f"{len(same)} files issued again by NOAA with the same content kept from earlier sets")
    log(f"set {tag}: {len(rows)} files, {total / 1e6:.0f} MB; new {len(new)} files, {added / 1e6:.0f} MB, "
        f"in {time.time() - started:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
