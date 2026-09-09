"""The set manifest: ``sets.txt`` parsing and manifest-driven ``split_ie`` dispatch.

The manifest exists so the saved set axis stops being read by component count. These
tests pin the two invariants that matter: on today's layouts the manifest must give
*exactly* what the count-based path gives (so shipping it changes nothing), and on the
layouts it is meant to enable -- continuum first, extra ladder sets -- it must still
resolve the three physics components correctly.
"""
import numpy as np
import pytest

from tumult.flux import ladder_sets, read_set_manifest, split_ie, write_set_manifest

HEADER = "# set  kind       element  da\n"


def write_manifest(folder, rows):
    (folder / 'sets.txt').write_text(HEADER + "".join(rows))
    return read_set_manifest(folder)


def fake_ie(ncomp, nlam=5, nang=3):
    """Distinct, non-degenerate components so a mislabelled slot cannot pass."""
    return np.arange(nlam * nang * ncomp, dtype=float).reshape(nlam, 1, 1, nang, ncomp)


def test_absent_manifest_is_none(tmp_path):
    assert read_set_manifest(tmp_path) is None


@pytest.mark.parametrize("rows, ncomp", [
    (["1 nlte - -\n", "2 continuum - -\n", "3 lte - -\n"], 3),
    (["1 lte - -\n", "2 continuum - -\n"], 2),
])
def test_matches_count_based_path_on_current_layouts(tmp_path, rows, ncomp):
    """Shipping the manifest must be a no-op for runs the old reader already handled."""
    sets = write_manifest(tmp_path, rows)
    ie = fake_ie(ncomp)
    by_manifest = split_ie(ie, False, sets)
    by_count = split_ie(ie, False)
    for a, b in zip(by_manifest[:3], by_count[:3]):
        assert np.array_equal(a, b)
    assert by_manifest[3] == by_count[3]


def test_lte_run_reports_lte_mode(tmp_path):
    sets = write_manifest(tmp_path, ["1 lte - -\n", "2 continuum - -\n"])
    ie = fake_ie(2)
    inten, inten_lte, _, mode = split_ie(ie, False, sets)
    assert mode == 'LTE'
    assert np.array_equal(inten, inten_lte)


def test_continuum_first_with_ladder_sets(tmp_path):
    """The layout the reorder and the ladder produce: position carries no meaning."""
    sets = write_manifest(tmp_path, [
        "1 continuum - -\n",
        "2 nlte - -\n",
        "3 lte - -\n",
        "4 nlte Fe -0.40\n",
        "5 nlte Fe -0.40\n",
        "5 nlte Ni +0.20\n",
    ])
    assert [s['kind'] for s in sets] == ['continuum', 'nlte', 'lte', 'nlte', 'nlte']
    # a multi-element set is several rows sharing one set number
    assert sets[4]['offsets'] == {'Fe': -0.40, 'Ni': 0.20}

    ie = fake_ie(5)
    inten, inten_lte, inten_cnt, mode = split_ie(ie, False, sets)
    assert np.array_equal(inten_cnt, ie[:, 0, 0, :, 0])
    assert np.array_equal(inten, ie[:, 0, 0, :, 1])
    assert np.array_equal(inten_lte, ie[:, 0, 0, :, 2])
    assert mode == 'NLTE'


def test_manifest_disagreeing_with_array_fails_loudly(tmp_path):
    sets = write_manifest(tmp_path, ["1 nlte - -\n", "2 continuum - -\n", "3 lte - -\n"])
    with pytest.raises(ValueError, match="manifest"):
        split_ie(fake_ie(2), False, sets)


def test_manifest_without_continuum_fails_loudly(tmp_path):
    sets = write_manifest(tmp_path, ["1 nlte - -\n", "2 lte - -\n"])
    with pytest.raises(ValueError, match="continuum"):
        split_ie(fake_ie(2), False, sets)


def test_base_set_indices_serves_other_arrays(tmp_path):
    """``i3`` shares the set axis, so it is unpacked by the same lookup."""
    from tumult.flux import base_set_indices

    sets = write_manifest(tmp_path, [
        "1 continuum - -\n", "2 nlte - -\n", "3 lte - -\n",
    ])
    assert base_set_indices(sets, 3) == (1, 2, 0, 'NLTE')

    sets = write_manifest(tmp_path, ["1 continuum - -\n", "2 lte - -\n"])
    assert base_set_indices(sets, 2) == (1, 1, 0, 'LTE')


def test_ladder_rows_and_writer_round_trip(tmp_path):
    sets = [
        {'kind': 'continuum', 'offsets': {}},
        {'kind': 'nlte', 'offsets': {}},
        {'kind': 'lte', 'offsets': {}},
        {'kind': 'nlte', 'offsets': {'Mg': -0.2, 'Si': -0.2}},
        {'kind': 'nlte', 'offsets': {'Mg': 0.2, 'Si': 0.2}},
    ]
    path = tmp_path / 'nested' / 'sets.txt'
    write_set_manifest(path, sets)
    assert read_set_manifest(path.parent) == sets
    assert ladder_sets(sets) == [
        (3, 'Mg', -0.2), (3, 'Si', -0.2),
        (4, 'Mg', 0.2), (4, 'Si', 0.2),
    ]


def test_manifest_rejects_conflicting_repeated_set(tmp_path):
    (tmp_path / 'sets.txt').write_text(
        HEADER + '1 continuum - -\n1 lte Fe +0.1\n')
    with pytest.raises(ValueError, match='conflicting'):
        read_set_manifest(tmp_path)
