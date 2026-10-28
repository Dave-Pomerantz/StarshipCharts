"""The depth grid file: one NOAA S-102 file converted for the devices.

Every cell keeps three values from NOAA's file: the depth, NOAA's stated uncertainty of that depth, and the
survey behind it (NOAA's QualityOfBathymetryCoverage). Coarser levels keep, for each block of cells, the
cell with the shallowest depth, and that cell's uncertainty and survey travel with it.

Layout, little-endian:
  header      magic b"STARGRID", u16 format (1), u16 EPSG of the grid's coordinates, u16 vertical datum
              (NOAA's code; 12 is mean lower low water), u16 block size, u16 level count, u16 reserved,
              f64 origin x, f64 origin y (the centre of cell 0,0), f64 spacing x, f64 spacing y,
              u32 columns, u32 rows (level 0), u32 survey table bytes
  surveys     UTF-8 text, tab separated, a header line then one line per survey; line 1 is survey index 1
  per level   u32 columns, u32 rows, then one entry per block, block rows south to north, block columns west
              to east: u64 offset, u32 length (0: the block holds no depth)
  blocks      zlib of: the block*block i16 depths (centimetres, positive down, -32768 empty), each row held as
              differences from the cell to its west (the first cell as is, wrapping at 16 bits), all low
              bytes then all high bytes; then block*block u16 uncertainty (decimetres, rounded up so it is
              never smaller than NOAA's, 65535 empty); then block*block u16 survey index (0 none)
Row 0 is the southern row; a level L cell covers level 0 cells [i*2^L, (i+1)*2^L) on each axis.
"""
import struct
import zlib

import numpy as np

MAGIC = b"STARGRID"
FORMAT = 1
BLOCK = 256
EMPTY_DEPTH = -32768
EMPTY_UNCERTAINTY = 65535
FILL = 1000000.0
HEADER = struct.Struct("<8sHHHHHHddddIII")
LEVEL = struct.Struct("<II")
ENTRY = struct.Struct("<QI")

SURVEY_FIELDS = [
    ("id", "id"),
    ("dateStart", "surveyDateRange.dateStart"),
    ("dateEnd", "surveyDateRange.dateEnd"),
    ("sourceSurveyID", "sourceSurveyID"),
    ("surveyAuthority", "surveyAuthority"),
    ("dataAssessment", "dataAssessment"),
    ("fullSeafloorCoverageAchieved", "fullSeafloorCoverageAchieved"),
    ("bathyCoverage", "bathyCoverage"),
    ("significantFeaturesDetected", "featuresDetected.significantFeaturesDetected"),
    ("leastDepthOfDetectedFeaturesMeasured", "featuresDetected.leastDepthOfDetectedFeaturesMeasured"),
    ("sizeOfFeaturesDetected", "featuresDetected.sizeOfFeaturesDetected"),
    ("featureSizeVar", "featureSizeVar"),
    ("horizontalUncertaintyFixed", "zoneOfConfidence.horizontalPositionUncertainty.uncertaintyFixed"),
    ("horizontalUncertaintyVariableFactor", "zoneOfConfidence.horizontalPositionUncertainty.uncertaintyVariableFactor"),
    ("typeOfBathymetricEstimationUncertainty", "typeOfBathymetricEstimationUncertainty"),
]


def _text(value):
    if isinstance(value, bytes):
        value = value.decode("utf-8", "replace")
    if isinstance(value, (float, np.floating)):
        value = repr(float(value))
    return str(value).replace("\t", " ").replace("\n", " ").replace("\r", " ")


def read_s102(path):
    """Reads one S-102 edition 3.0 file: (attributes, depth cm i16, uncertainty dm u16, survey index u16, survey lines)."""
    import h5py

    with h5py.File(path, "r") as f:
        cov = f["BathymetryCoverage/BathymetryCoverage.01"]
        values = cov["Group_001/values"][:]
        quality = f["QualityOfBathymetryCoverage/QualityOfBathymetryCoverage.01/Group_001/values"][:]
        table = f["QualityOfBathymetryCoverage/featureAttributeTable"][:]
        attrs = {
            "epsg": int(f.attrs["horizontalCRS"]),
            "datum": int(f.attrs["verticalDatum"]),
            "ox": float(cov.attrs["gridOriginLongitude"]),
            "oy": float(cov.attrs["gridOriginLatitude"]),
            "sx": float(cov.attrs["gridSpacingLongitudinal"]),
            "sy": float(cov.attrs["gridSpacingLatitudinal"]),
        }
    depth = values["depth"]
    unc = values["uncertainty"]
    valid = depth != FILL
    depth_cm = np.where(valid, np.clip(np.round(depth.astype(np.float64) * 100), -32767, 32767), EMPTY_DEPTH).astype("<i2")
    unc_ok = valid & (unc != FILL)
    # Decimetres, rounded up: never smaller than NOAA's figure, and a third of the size of centimetres.
    unc_dm = np.where(unc_ok, np.clip(np.ceil(np.round(unc.astype(np.float64) * 100) / 10), 0, 65534), EMPTY_UNCERTAINTY).astype("<u2")
    # NOAA's survey ids are large numbers; the file keeps an index into its own table instead (1 is the table's
    # first line), so a cell holds two bytes.
    ids = [int(r["id"]) for r in table]
    if len(ids) > 65534:
        raise ValueError(f"{path}: {len(ids)} surveys do not fit the survey index")
    lookup = np.zeros(max(ids, default=0) + 1, dtype="<u2")
    for i, sid in enumerate(ids):
        lookup[sid] = i + 1
    q = quality.astype(np.int64)
    known = (q > 0) & (q < lookup.size)
    survey = np.where(valid & known, lookup[np.where(known, q, 0)], 0).astype("<u2")
    names = table.dtype.names
    lines = ["\t".join(k for k, _ in SURVEY_FIELDS)]
    for r in table:
        lines.append("\t".join(_text(r[src]) if src in names else "" for _, src in SURVEY_FIELDS))
    return attrs, depth_cm, unc_dm, survey, lines


def _stack(a):
    return np.stack([a[0::2, 0::2], a[0::2, 1::2], a[1::2, 0::2], a[1::2, 1::2]])


def coarser(depth, unc, survey):
    """One level up: in each 2 by 2 group the shallowest depth, with its own uncertainty and survey."""
    rows, cols = depth.shape
    r2, c2 = (rows + 1) // 2, (cols + 1) // 2

    def pad(a, v):
        return np.pad(a, ((0, r2 * 2 - rows), (0, c2 * 2 - cols)), constant_values=v)

    big = 1 << 20
    d = pad(depth, EMPTY_DEPTH).astype(np.int32)
    d = np.where(d == EMPTY_DEPTH, big, d)
    ds, us, ss = _stack(d), _stack(pad(unc, EMPTY_UNCERTAINTY)), _stack(pad(survey, 0))
    pick = np.argmin(ds, axis=0)[None]
    dd = np.take_along_axis(ds, pick, axis=0)[0]
    out_d = np.where(dd == big, EMPTY_DEPTH, dd).astype("<i2")
    return out_d, np.take_along_axis(us, pick, axis=0)[0].astype("<u2"), np.take_along_axis(ss, pick, axis=0)[0].astype("<u2")


def encode_depth(a):
    """Each row as differences from the cell to its west, wrapping at 16 bits; low bytes, then high bytes."""
    v = a.view(np.uint16).astype(np.int64)
    diff = v.copy()
    diff[:, 1:] = (v[:, 1:] - v[:, :-1]) & 0xFFFF
    planes = diff.astype("<u2").view(np.uint8).reshape(a.shape + (2,))
    return planes[..., 0].tobytes() + planes[..., 1].tobytes()


def decode_depth(raw, block):
    k = block * block
    low = np.frombuffer(raw, np.uint8, k, 0).astype(np.int64)
    high = np.frombuffer(raw, np.uint8, k, k).astype(np.int64)
    diff = (low | (high << 8)).reshape(block, block)
    return (np.cumsum(diff, axis=1) & 0xFFFF).astype("<u2").view("<i2")


def _full(a, block, value, dtype):
    return np.pad(a, ((0, block - a.shape[0]), (0, block - a.shape[1])), constant_values=value).astype(dtype)


def write_grid(path, attrs, depth, unc, survey, lines, block=BLOCK):
    levels = [(depth, unc, survey)]
    while max(levels[-1][0].shape) > block:
        levels.append(coarser(*levels[-1]))
    table = ("\n".join(lines) + "\n").encode("utf-8")
    rows, cols = depth.shape
    index = []
    for d, u, s in levels:
        lr, lc = d.shape
        br, bc = -(-lr // block), -(-lc // block)
        entries = []
        for by in range(br):
            for bx in range(bc):
                sl = (slice(by * block, (by + 1) * block), slice(bx * block, (bx + 1) * block))
                bd = d[sl]
                if not (bd != EMPTY_DEPTH).any():
                    entries.append(None)
                    continue
                raw = (encode_depth(_full(bd, block, EMPTY_DEPTH, "<i2"))
                       + _full(u[sl], block, EMPTY_UNCERTAINTY, "<u2").tobytes()
                       + _full(s[sl], block, 0, "<u2").tobytes())
                entries.append(zlib.compress(raw, 6))
        index.append((lc, lr, entries))
    offset = HEADER.size + len(table) + sum(LEVEL.size + ENTRY.size * len(e) for _, _, e in index)
    with open(path, "wb") as f:
        f.write(HEADER.pack(MAGIC, FORMAT, attrs["epsg"], attrs["datum"], block, len(levels), 0,
                            attrs["ox"], attrs["oy"], attrs["sx"], attrs["sy"], cols, rows, len(table)))
        f.write(table)
        for lc, lr, entries in index:
            f.write(LEVEL.pack(lc, lr))
            for blob in entries:
                if blob is None:
                    f.write(ENTRY.pack(0, 0))
                else:
                    f.write(ENTRY.pack(offset, len(blob)))
                    offset += len(blob)
        for _, _, entries in index:
            for blob in entries:
                if blob is not None:
                    f.write(blob)


def convert(h5_path, grid_path):
    attrs, depth, unc, survey, lines = read_s102(h5_path)
    write_grid(grid_path, attrs, depth, unc, survey, lines)
    return attrs


def read_grid(path):
    """Reads a whole grid file back, for checks: the header and, per level, (depth, uncertainty, survey) arrays."""
    with open(path, "rb") as f:
        data = f.read()
    h = HEADER.unpack_from(data, 0)
    if h[0] != MAGIC:
        raise ValueError("not a depth grid file")
    keys = ["magic", "format", "epsg", "datum", "block", "levels", "reserved", "ox", "oy", "sx", "sy", "cols", "rows", "table"]
    head = dict(zip(keys, h))
    pos = HEADER.size
    head["surveys"] = data[pos:pos + head["table"]].decode("utf-8").splitlines()
    pos += head["table"]
    block = head["block"]
    k = block * block
    levels = []
    for _ in range(head["levels"]):
        lc, lr = LEVEL.unpack_from(data, pos)
        pos += LEVEL.size
        br, bc = -(-lr // block), -(-lc // block)
        d = np.full((br * block, bc * block), EMPTY_DEPTH, "<i2")
        u = np.full(d.shape, EMPTY_UNCERTAINTY, "<u2")
        s = np.zeros(d.shape, "<u2")
        for by in range(br):
            for bx in range(bc):
                off, n = ENTRY.unpack_from(data, pos)
                pos += ENTRY.size
                if n == 0:
                    continue
                raw = zlib.decompress(data[off:off + n])
                sl = (slice(by * block, (by + 1) * block), slice(bx * block, (bx + 1) * block))
                d[sl] = decode_depth(raw[:2 * k], block)
                u[sl] = np.frombuffer(raw, "<u2", k, 2 * k).reshape(block, block)
                s[sl] = np.frombuffer(raw, "<u2", k, 4 * k).reshape(block, block)
        levels.append((d[:lr, :lc], u[:lr, :lc], s[:lr, :lc]))
    return head, levels
