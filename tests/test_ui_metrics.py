"""Tests for metric hierarchy selection used by the Streamlit UI."""

from gabes.ui_metrics import partition_metrics, split_metric_value


def test_metric_value_splits_number_from_single_token_unit():
    assert split_metric_value("-59.89 dBm") == ("-59.89", "dBm")
    assert split_metric_value("+2.264e-4 /MHz") == ("+2.264e-4", "/MHz")
    assert split_metric_value("800.0 µA") == ("800.0", "µA")
    assert split_metric_value(" −0.25 mrad/µT ") == ("−0.25", "mrad/µT")


def test_metric_value_keeps_unitless_numbers_and_status_text_intact():
    values = (
        15.2,
        "invalid transition",
        "EIA peak",
        "calibrated · non-predictive",
        "404 ERROR state",
        "—",
        "<script>alert('metric')</script>",
    )

    for value in values:
        assert split_metric_value(value) == (str(value).strip(), None)


def test_metric_value_never_splits_explicit_status_text():
    assert split_metric_value("2 peaks", kind="status") == ("2 peaks", None)
    assert split_metric_value("404 ERROR", kind="status") == ("404 ERROR", None)


def test_explicit_heroes_take_priority_and_missing_slot_is_filled_in_order():
    first = {"label": "first", "value": 1}
    explicit = {"label": "explicit", "value": 2, "tier": "hero"}
    last = {"label": "last", "value": 3}
    metrics = [first, explicit, last]

    heroes, ribbon = partition_metrics(metrics)

    assert heroes == [explicit, first]
    assert ribbon == [last]
    assert heroes[0] is explicit
    assert heroes[1] is first


def test_first_metrics_become_heroes_when_none_are_tagged():
    metrics = [
        {"label": "one", "value": 1},
        {"label": "two", "value": 2},
        {"label": "three", "value": 3},
    ]

    heroes, ribbon = partition_metrics(metrics)

    assert heroes == metrics[:2]
    assert ribbon == metrics[2:]


def test_explicit_hero_overflow_returns_to_ribbon_in_original_order():
    hero_1 = {"label": "hero 1", "value": 1, "tier": "hero"}
    ordinary = {"label": "ordinary", "value": 2}
    hero_2 = {"label": "hero 2", "value": 3, "tier": "hero"}
    hero_3 = {"label": "hero 3", "value": 4, "tier": "hero"}
    tail = {"label": "tail", "value": 5}

    heroes, ribbon = partition_metrics(
        [hero_1, ordinary, hero_2, hero_3, tail], hero_count=2
    )

    assert heroes == [hero_1, hero_2]
    assert ribbon == [ordinary, hero_3, tail]


def test_empty_singleton_and_zero_hero_count_are_deterministic():
    assert partition_metrics([]) == ([], [])

    only = {"label": "only", "value": 1}
    heroes, ribbon = partition_metrics([only])
    assert heroes == [only]
    assert ribbon == []

    heroes, ribbon = partition_metrics([only], hero_count=0)
    assert heroes == []
    assert ribbon == [only]


def test_partition_does_not_mutate_metrics_or_input_list():
    explicit = {"label": "explicit", "value": 1, "tier": "hero"}
    ordinary = {"label": "ordinary", "value": 2, "meta": {"unit": "MHz"}}
    metrics = [ordinary, explicit]
    original_items = list(metrics)
    original_dicts = [dict(metric) for metric in metrics]

    heroes, ribbon = partition_metrics(metrics, hero_count=1)

    assert metrics == original_items
    assert metrics[0] is ordinary
    assert metrics[1] is explicit
    assert metrics == original_dicts
    assert heroes[0] is explicit
    assert ribbon[0] is ordinary
