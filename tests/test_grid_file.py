"""Round trip of the depth grid file on a small made-up S-102 file: every cell's depth, uncertainty and survey
come back, and each coarser level holds the shallowest cell of its group with that cell's uncertainty and survey."""
import os
import sys
import tempfile

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "chartset"))
import grid_file as g  # noqa: E402


def make_s102(path, rows, cols, seed=1):
    rng = np.random.default_rng(seed)
    depth = rng.uniform(-2, 40, (rows, cols)).astype("<f4")
    unc = rng.uniform(0.1, 3, (rows, cols)).astype("<f4")
    depth[0:5, 0:7] = g.FILL
    unc[0:5, 0:7] = g.FILL
    ids = np.array([0, 108050, 7], dtype="<u4")
    quality = rng.choice(ids[1:], (rows, cols)).astype("<u4")
    quality[0:5, 0:7] = 0
    with h5py.File(path, "w") as f:
        f.attrs["horizontalCRS"] = 32619
        f.attrs["verticalDatum"] = 12
        cov = f.create_group("BathymetryCoverage/BathymetryCoverage.01")
        for k, v in [("gridOriginLongitude", 351511.5), ("gridOriginLatitude", 4676200.5), ("gridSpacingLongitudinal", 4.0), ("gridSpacingLatitudinal", 4.0)]:
            cov.attrs[k] = v
        vals = np.zeros((rows, cols), dtype=[("depth", "<f4"), ("uncertainty", "<f4")])
        vals["depth"], vals["uncertainty"] = depth, unc
        cov.create_dataset("Group_001/values", data=vals)
        f.create_dataset("QualityOfBathymetryCoverage/QualityOfBathymetryCoverage.01/Group_001/values", data=quality)
        table = np.zeros(2, dtype=[("id", "<u4"), ("dataAssessment", "u1"), ("surveyDateRange.dateStart", h5py.string_dtype()),
                                   ("sourceSurveyID", h5py.string_dtype())])
        table["id"] = [108050, 7]
        table["dataAssessment"] = [1, 3]
        table["surveyDateRange.dateStart"] = ["2024-08-23", "1807-02-10"]
        table["sourceSurveyID"] = ["H13997", "Generalization"]
        f.create_dataset("QualityOfBathymetryCoverage/featureAttributeTable", data=table)
    return depth, unc, quality


def main():
    with tempfile.TemporaryDirectory() as tmp:
        src, out = os.path.join(tmp, "t.h5"), os.path.join(tmp, "t.grid")
        depth, unc, quality = make_s102(src, 600, 530)
        g.convert(src, out)
        head, levels = g.read_grid(out)
        assert head["epsg"] == 32619 and head["datum"] == 12 and head["cols"] == 530 and head["rows"] == 600
        assert len(levels) == 3, len(levels)
        d, u, s = levels[0]
        valid = depth != g.FILL
        assert (d[~valid] == g.EMPTY_DEPTH).all() and (u[~valid] == g.EMPTY_UNCERTAINTY).all() and (s[~valid] == 0).all()
        assert np.array_equal(d[valid], np.round(depth[valid].astype(np.float64) * 100).astype(np.int16))
        assert (u[valid] * 10 >= np.round(unc[valid].astype(np.float64) * 100)).all(), "uncertainty never rounded down"
        assert (u[valid] * 10 - np.round(unc[valid].astype(np.float64) * 100) < 10).all()
        assert np.array_equal(s[valid], np.where(quality[valid] == 108050, 1, 2))
        assert head["surveys"][1].startswith("108050\t2024-08-23") and head["surveys"][2].startswith("7\t1807-02-10")
        for level in range(1, len(levels)):
            expected = g.coarser(*levels[level - 1])
            assert all(np.array_equal(a, b) for a, b in zip(levels[level], expected))
            fine, coarse = levels[level - 1][0].astype(np.int64), levels[level][0].astype(np.int64)
            fine = np.where(fine == g.EMPTY_DEPTH, 1 << 20, fine)
            for r, c in [(0, 0), (7, 9), (coarse.shape[0] - 1, coarse.shape[1] - 1)]:
                group = fine[2 * r:2 * r + 2, 2 * c:2 * c + 2]
                want = group.min()
                assert coarse[r, c] == (g.EMPTY_DEPTH if want == 1 << 20 else want)
    print("grid file round trip: ok")


if __name__ == "__main__":
    main()
