"""Regression tests for batched scientific updates."""
from types import SimpleNamespace


def test_local_continuum_batch_aggregates_once(monkeypatch):
    from tumult import abundance as ab

    lines = [SimpleNamespace(value=0), SimpleNamespace(value=0)]

    class Result:
        def __init__(self):
            self.lines = lines
            self.aggregations = 0

        def reaggregate(self):
            self.aggregations += 1

    def edit(result, line, target, **kwargs):
        return SimpleNamespace(value=kwargs["value"])

    monkeypatch.setattr(ab, "_local_continuum_edit", edit)
    result = Result()
    failed = ab.set_local_continua(
        result, [(lines[0], {"value": 1}), (lines[1], {"value": 2})], object())

    assert not failed
    assert [line.value for line in lines] == [1, 2]
    assert result.aggregations == 1
