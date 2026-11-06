"""Native comparison rows, bounded row cache and cancellable sequential sweeps.

No physics or interpolation lives here. Every new row is obtained from the
scheme's own comparison adapter and headless readout.
"""
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import hashlib
import json
import math
from threading import RLock
from time import monotonic

from gabes.model_comparison import comparison_values

_EXECUTOR = ThreadPoolExecutor(max_workers=1, thread_name_prefix="gabes-spectrum")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def workspace_key(scheme, params, varied_key, readout_version):
    return digest([scheme.name, scheme.cache_version, readout_version, varied_key,
                   {k: v for k, v in params.items() if k != varied_key}])


def row_value(row):
    return row.get("parameter", {}).get("value")


def reusable_dataset(entry, scheme, params, readout_version):
    """Selecting an existing row must not recenter or repeat the comparison."""
    if not entry:
        return None
    data = entry["dataset"]
    key = data.get("varied", {}).get("key")
    if not key or key not in params:
        return None
    if entry["identity"] != workspace_key(scheme, params, key, readout_version):
        entry["runner"].cancel()
        return None
    if any(row_value(row) == params[key] for row in data["series"]):
        return data
    entry["runner"].cancel()
    return None


def reusable_solution(entry, scheme, params, readout_version):
    if not entry:
        return None
    key = entry["dataset"].get("varied", {}).get("key")
    if not key or key not in params or entry["identity"] != workspace_key(scheme, params, key, readout_version):
        entry["runner"].cancel()
        return None
    return entry["runner"].solution(entry["identity"], params[key])


def cancel_inactive_sweeps(session, active_scheme):
    for key in list(session):
        if key.startswith("spectral-workspace-") and key != "spectral-workspace-"+active_scheme:
            session[key]["runner"].cancel()


def extract_row(dataset, index):
    """Preserve native coordinates, gaps, scope and native summary coordinates."""
    source = dataset["series"][index]
    row = deepcopy(source)
    parameter = row.get("parameter", {})
    row["id"] = "condition-" + digest([parameter.get("key"), parameter.get("value")])[:16]
    panels = []
    for panel in dataset["panels"]:
        curves = [dict(c, series_id=row["id"]) for c in panel["curves"]
                  if c["series_id"] == source["id"]]
        if curves:
            panels.append({**deepcopy(panel), "curves": curves})
    result = {"series": row, "panels": panels, "axis": deepcopy(dataset["axis"])}
    summary = dataset.get("summary")
    if summary and all(len(s["x"]) == len(dataset["series"]) and
                       all(len(y) == len(s["x"]) for y in s["observables"].values())
                       for s in summary["series"]):
        result["summary"] = deepcopy(summary)
        for curve in result["summary"]["series"]:
            curve["x"] = [curve["x"][index]]
            curve["observables"] = {k: [v[index]] for k, v in curve["observables"].items()}
            for field in ("sample_status", "sample_reasons"):
                if field in curve:
                    curve[field] = [curve[field][index]]
            if "observable_status" in curve:
                curve["observable_status"] = {k: [v[index]] for k, v in curve["observable_status"].items()}
            curve["condition_ids"] = [row["id"]]
    return result


def resolve_row(scheme, params, key, value, resolve_raw):
    condition = {**params, key: value}
    raw = resolve_raw(condition)
    view = scheme.headless_observables(raw, condition)
    with comparison_values(key, [value]):
        data = scheme.model_comparison(condition, resolve_raw=resolve_raw,
                                       current_raw=raw, current_view=view)
    if data.get("schema") != "gabes-spectral-comparison-v1":
        raise ValueError("This comparison does not expose native spectral rows.")
    matches = [i for i, row in enumerate(data["series"])
               if row.get("parameter", {}).get("key") == key and row_value(row) == value]
    if not matches:
        raise ValueError("The scheme did not return the requested comparison condition.")
    result = extract_row(data, matches[0])
    result.update(_native_raw=raw, _native_view=view)
    return result


def scan_values(low, high, count, current, spec):
    """Explicit uniform samples plus the exact condition at scan start."""
    numbers = (low, high, count, current)
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or
           not math.isfinite(x) for x in numbers):
        raise ValueError("Enter finite scan limits and a whole sample count.")
    if int(count) != count or not 2 <= count <= 41:
        raise ValueError("Choose 2–41 uniform condition samples.")
    if low >= high:
        raise ValueError("The upper scan limit must exceed the lower limit.")
    if ((spec.vmin is not None and low < spec.vmin) or
            (spec.vmax is not None and high > spec.vmax)):
        raise ValueError("Keep the condition scan inside the parameter's supported range.")
    values = [low + (high-low)*i/(int(count)-1) for i in range(int(count))]
    if isinstance(spec.default, int) and not isinstance(spec.default, bool):
        values = [int(round(v)) for v in values]
    values += [current]
    return sorted(set(values))


def assemble_dataset(template, snapshot):
    data = deepcopy(template)
    rows = snapshot["rows"]
    if not rows:
        return data
    data["series"] = [deepcopy(payload["series"]) for _, payload in rows]
    panels = OrderedDict()
    for _, payload in rows:
        if payload["axis"] != template["axis"]:
            raise ValueError("A scan row changed the physical axis or its reference.")
        for panel in payload["panels"]:
            target = panels.setdefault(panel["key"], {**deepcopy(panel), "curves": []})
            target["curves"].extend(deepcopy(panel["curves"]))
    data["panels"] = list(panels.values())
    summaries = [payload.get("summary") for _, payload in rows]
    if all(summaries):
        summary = deepcopy(summaries[0])
        for curve in summary["series"]:
            others = [next((s for s in part["series"] if s["id"] == curve["id"]), None)
                      for part in summaries]
            if any(s is None for s in others):
                raise ValueError("Summary rows changed their observable contract.")
            curve["x"] = [x for s in others for x in s["x"]]
            curve["observables"] = {k: [v for s in others for v in s["observables"][k]]
                                     for k in curve["observables"]}
            curve["condition_ids"] = [v for s in others for v in s["condition_ids"]]
            for field in ("sample_status", "sample_reasons"):
                if all(field in s for s in others):
                    curve[field] = [v for s in others for v in s[field]]
            if all("observable_status" in s for s in others):
                curve["observable_status"] = {k: [v for s in others for v in s["observable_status"][k]]
                                              for k in curve["observable_status"]}
        data["summary"] = summary
    else:
        data.pop("summary", None)
    data["varied"]["values"] = [row_value(row) for row in data["series"]]
    data["sweep"] = {k: deepcopy(v) for k, v in snapshot.items() if k != "rows"}
    return data


class RowSweep:
    """One in-flight row at a time; cancellation prevents subsequent rows.

    Completed rows remain usable after cancellation or failure. A stale job
    cannot publish into a replacement job. Cache entries include all fixed
    inputs and display/readout versions in the supplied fingerprint.
    """
    def __init__(self, max_entries=64, lease_seconds=120):
        self._lock = RLock()
        self._cache = OrderedDict()
        self._maximum = max_entries
        self._lease = lease_seconds
        self._last_seen = monotonic()
        self._generation = 0
        self._job = None

    def _remember(self, key, payload):
        self._cache[key] = deepcopy(payload)
        self._cache.move_to_end(key)
        while len(self._cache) > self._maximum:
            self._cache.popitem(last=False)

    def seed(self, fingerprint, values_and_rows):
        with self._lock:
            for value, payload in values_and_rows:
                self._remember((fingerprint, value), payload)

    def solution(self, fingerprint, value):
        """Current plotted rows retain solves even if the global cache evicts them."""
        with self._lock:
            job = self._job
            payload = (job["rows"].get(value) if job and job["fingerprint"] == fingerprint else None)
            payload = payload or self._cache.get((fingerprint, value))
            if payload is None or "_native_raw" not in payload:
                return None
            return deepcopy((payload["_native_raw"], payload.get("_native_view")))

    def start(self, fingerprint, values, resolver, *, request=None):
        with self._lock:
            self._generation += 1
            generation = self._generation
            values = tuple(dict.fromkeys(values))
            rows = {value: deepcopy(self._cache[(fingerprint, value)])
                    for value in values if (fingerprint, value) in self._cache}
            self._last_seen = monotonic()
            self._resolver = resolver
            self._job = dict(generation=generation, fingerprint=fingerprint,
                             requested_values=values, rows=rows, status="running",
                             error=None, active_value=None, cached=len(rows), request=deepcopy(request))
        _EXECUTOR.submit(self._run, generation, resolver)

    def _run(self, generation, resolver):
        while True:
            with self._lock:
                job = self._job
                if job is None or generation != job["generation"] or job["status"] != "running":
                    return
                if monotonic()-self._last_seen > self._lease:
                    job["status"] = "cancelled"
                    return
                pending = [v for v in job["requested_values"] if v not in job["rows"]]
                if not pending:
                    job["status"] = "complete"
                    job["active_value"] = None
                    return
                value = pending[0]
                fingerprint = job["fingerprint"]
                cached = self._cache.get((fingerprint, value))
                if cached is not None:
                    job["rows"][value] = deepcopy(cached)
                    job["cached"] += 1
                    continue
                job["active_value"] = value
            try:
                payload = resolver(value)
            except Exception as error:
                with self._lock:
                    if self._job is job and job["status"] == "running":
                        job.update(status="failed", error=str(error), active_value=None)
                return
            with self._lock:
                self._remember((fingerprint, value), payload)
                if self._job is not job or job["status"] != "running":
                    return
                job["rows"][value] = payload

    def cancel(self):
        with self._lock:
            if self._job is not None and self._job["status"] == "running":
                self._job.update(status="cancelled", active_value=None)

    def resume(self):
        with self._lock:
            job = self._job
            if job is None or job["status"] not in {"cancelled", "failed"}:
                return
            self.seed(job["fingerprint"], job["rows"].items())
            self.start(job["fingerprint"], job["requested_values"], self._resolver, request=job["request"])

    def snapshot(self):
        with self._lock:
            self._last_seen = monotonic()
            if self._job is None:
                return None
            job = self._job
            return dict(generation=job["generation"], status=job["status"],
                        requested_values=list(job["requested_values"]),
                        completed=len(job["rows"]), total=len(job["requested_values"]),
                        cached=job["cached"], active_value=job["active_value"], error=job["error"],
                        request=deepcopy(job["request"]),
                        rows=[(v, deepcopy({k: val for k, val in job["rows"][v].items()
                                            if not k.startswith("_native_")})) for v in job["requested_values"]
                              if v in job["rows"]])
