"""Exact Python-text provenance contracts for new campaign runs."""

import builtins
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import runpy

import pytest

from analysis.grand_challenge import source_provenance as provenance
from analysis.grand_challenge.source_provenance import (
    capture_sources, require_sources, source_identity,
)


SOURCE = ('# coding: utf-8\n'
          '# Preserve comments and coefficients.\n'
          'COEFFICIENT = 0.125\n'
          'LABEL = "Rb α"\n'
          'TEXT = """first\nsecond"""\n'
          'ESCAPED = "literal \\r\\n"\n')
PATH = 'physics/model.py'


def write_source(root, name=PATH, raw=None):
    path = root.joinpath(*name.split('/'))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(SOURCE.encode('utf-8') if raw is None else raw)
    return path


@pytest.fixture
def captured(tmp_path):
    path = write_source(tmp_path)
    return tmp_path, path, capture_sources(tmp_path, [PATH])


def test_relocation_and_lf_crlf_return_fresh_raw_audit(tmp_path):
    original, relocated = tmp_path/'original', tmp_path/'relocated'
    lf = SOURCE.encode('utf-8')
    crlf = lf.replace(b'\n', b'\r\n')
    write_source(original, raw=lf)
    write_source(relocated, raw=crlf)
    manifest = capture_sources(original, [PATH])
    saved = deepcopy(manifest)
    fresh = require_sources(relocated, json.loads(json.dumps(manifest)))
    assert manifest == saved
    assert source_identity(manifest) == source_identity(fresh)
    assert fresh == capture_sources(relocated, [PATH])
    assert fresh['manifest_sha256'] != manifest['manifest_sha256']
    old, new = manifest['files'][0], fresh['files'][0]
    assert old['path'] == new['path'] == PATH
    assert old['encoding'] == new['encoding'] == 'utf-8'
    assert old['raw_sha256'] == hashlib.sha256(lf).hexdigest()
    assert new['raw_sha256'] == hashlib.sha256(crlf).hexdigest()
    assert old['raw_sha256'] != new['raw_sha256']
    assert old['raw_bytes'] == len(lf)
    assert new['raw_bytes'] == len(crlf)
    assert old['normalized_sha256'] == new['normalized_sha256'] == hashlib.sha256(lf).hexdigest()
    assert old['normalized_bytes'] == new['normalized_bytes'] == len(lf)
    assert str(original) not in json.dumps(manifest)
    assert str(relocated) not in json.dumps(fresh)


@pytest.mark.parametrize('ending', [b'\r', b'\r\n', b'\n'])
def test_universal_newlines_include_multiline_strings(captured, ending):
    root, path, manifest = captured
    path.write_bytes(SOURCE.encode('utf-8').replace(b'\n', ending))
    fresh = require_sources(root, manifest)
    assert source_identity(fresh) == source_identity(manifest)


@pytest.mark.parametrize('ending', [b'\r', b'\r\n', b'\n'])
def test_encoding_cookie_detection_uses_only_first_two_physical_lines(tmp_path, ending):
    raw = b'# header\n# ordinary comment\nLABEL = "coding: latin-1"\n'
    write_source(tmp_path, raw=raw)
    manifest = capture_sources(tmp_path, [PATH])
    write_source(tmp_path, raw=raw.replace(b'\n', ending))
    assert source_identity(require_sources(tmp_path, manifest)) == source_identity(manifest)


def test_mixed_newlines_normalize_without_splitting_other_separators(tmp_path):
    raw = b'# header\r\nTEXT = "first\rsecond\nthird"\r'
    # Use triple quotes so all physical line endings are valid inside the string.
    raw = raw.replace(b'"', b'"""')
    write_source(tmp_path, raw=raw)
    manifest = capture_sources(tmp_path, [PATH])
    canonical = b'# header\nTEXT = """first\nsecond\nthird"""\n'
    assert manifest['files'][0]['normalized_sha256'] == hashlib.sha256(canonical).hexdigest()
    write_source(tmp_path, raw=canonical)
    assert source_identity(require_sources(tmp_path, manifest)) == source_identity(manifest)


@pytest.mark.parametrize('changed', [
    SOURCE.replace('0.125', '0.126'),
    SOURCE.replace('0.125', '1/8'),  # Numerically equal still changes source.
    SOURCE.replace('Rb α', 'Rb β'),
    SOURCE.replace('first\nsecond', 'first\nthird'),
    SOURCE.replace('literal \\r\\n', 'literal \\n'),
    SOURCE.replace('"Rb α"', "'Rb α'"),
    SOURCE.replace('Preserve comments', 'Keep comments'),
    SOURCE.replace('COEFFICIENT =', 'COEFFICIENT  ='),
    SOURCE.replace('0.125\n', '0.125 \n'),
    SOURCE + '\n',
    SOURCE.rstrip('\n'),
    SOURCE.replace('Rb α', 'Rb α\u2028'),
])
def test_any_source_edit_except_universal_newlines_changes_identity(captured, changed):
    root, path, manifest = captured
    path.write_bytes(changed.encode('utf-8'))
    assert source_identity(capture_sources(root, [PATH])) != source_identity(manifest)
    with pytest.raises(ValueError, match='source changed: physics/model.py'):
        require_sources(root, manifest)


def test_deletion_is_rejected(captured):
    root, path, manifest = captured
    identity = source_identity(manifest)
    path.unlink()
    assert source_identity(manifest) == identity  # Manifest validation is offline.
    with pytest.raises(ValueError, match='cannot read source file'):
        require_sources(root, manifest)


def test_only_explicit_paths_are_captured_in_stable_order(tmp_path):
    write_source(tmp_path, 'z.py')
    write_source(tmp_path, 'a.py')
    unlisted = write_source(tmp_path, 'unlisted.py', b'this is invalid Python !!!')
    manifest = capture_sources(tmp_path, ['z.py', 'a.py'])
    assert [entry['path'] for entry in manifest['files']] == ['a.py', 'z.py']
    assert capture_sources(tmp_path, iter(['a.py', 'z.py'])) == manifest
    unlisted.unlink()
    assert require_sources(tmp_path, manifest) == manifest
    assert source_identity(capture_sources(tmp_path, ['a.py'])) != source_identity(manifest)
    assert source_identity(capture_sources(tmp_path, ['a.py'])) != source_identity(
        capture_sources(tmp_path, ['z.py']))


@pytest.mark.parametrize('paths', [
    [], None, PATH, PATH.encode(), {PATH: 'legacy-sha'},
    ['missing.py'], [PATH, PATH], [PATH, 'PHYSICS/model.py'],
    ['../outside.py'], ['physics/../model.py'], ['/absolute.py'],
    ['C:/absolute.py'], ['C:relative.py'], ['//host/share/source.py'],
    ['physics\\model.py'], ['./physics/model.py'], ['physics//model.py'],
    ['physics/./model.py'], ['physics/model.py/'], ['physics/model.txt'],
    ['physics/model.PY'], ['physics/*.py'], ['physics/model.py:stream.py'],
    ['NUL.py'], ['CON/model.py'], ['trailing /model.py'], ['trailing./model.py'],
    ['bad\x00.py'], ['bad\n.py'], [Path(PATH)], [42],
])
def test_missing_unsafe_non_python_or_duplicate_paths_are_rejected(captured, paths):
    root, _, _ = captured
    with pytest.raises(ValueError):
        capture_sources(root, paths)


def test_missing_root_and_directories_are_rejected(tmp_path):
    with pytest.raises(ValueError, match='root'):
        capture_sources(tmp_path/'missing', ['model.py'])
    folder = tmp_path/'folder.py'
    folder.mkdir()
    with pytest.raises(ValueError, match='regular file'):
        capture_sources(tmp_path, ['folder.py'])
    file = write_source(tmp_path)
    with pytest.raises(ValueError, match='root'):
        capture_sources(file, [PATH])


def make_symlink(link, target, *, directory=False):
    try:
        link.symlink_to(target, target_is_directory=directory)
    except (OSError, NotImplementedError) as exc:
        pytest.skip(f'filesystem does not permit symlinks: {exc}')


@pytest.mark.parametrize('directory', [False, True])
def test_symlink_escape_is_rejected(tmp_path, directory):
    root, outside = tmp_path/'repo', tmp_path/'outside'
    root.mkdir()
    target = write_source(outside, 'outside.py')
    link = root/('linked' if directory else 'linked.py')
    make_symlink(link, outside if directory else target, directory=directory)
    with pytest.raises(ValueError, match='escapes root'):
        capture_sources(root, ['linked/outside.py' if directory else 'linked.py'])


def test_internal_symlink_is_allowed_but_alias_duplicates_are_rejected(captured):
    root, path, _ = captured
    make_symlink(root/'alias.py', path)
    manifest = capture_sources(root, ['alias.py'])
    assert require_sources(root, manifest) == manifest
    with pytest.raises(ValueError, match='duplicate'):
        capture_sources(root, ['alias.py', PATH])


def test_hardlink_alias_duplicates_are_rejected(captured):
    root, path, _ = captured
    try:
        (root/'alias.py').hardlink_to(path)
    except (OSError, NotImplementedError) as exc:
        pytest.skip(f'filesystem does not permit hardlinks: {exc}')
    with pytest.raises(ValueError, match='duplicate'):
        capture_sources(root, ['alias.py', PATH])


@pytest.mark.parametrize('raw', [
    b'def broken(:\n', b'return 1\n', b'value = 1\x00\n',
    b'label = "\xff"\n', b'# coding: unknown-codec\nx = 1\n',
    b'\xef\xbb\xbf# coding: latin-1\nx = 1\n',
])
def test_requires_decodable_compilable_python(tmp_path, raw):
    write_source(tmp_path, raw=raw)
    with pytest.raises(ValueError, match='invalid Python source or encoding'):
        capture_sources(tmp_path, [PATH])


@pytest.mark.parametrize('encoding', ['latin-1', 'utf-8-sig'])
def test_encoding_detection_and_utf8_normalized_byte_counts(tmp_path, encoding):
    text = '# coding: ' + ('latin-1' if encoding == 'latin-1' else 'utf-8') + '\nTEXT = "café"\n'
    raw = text.replace('\n', '\r\n').encode(encoding)
    write_source(tmp_path, raw=raw)
    manifest = capture_sources(tmp_path, [PATH])
    entry = manifest['files'][0]
    assert entry['encoding'] == ('iso-8859-1' if encoding == 'latin-1' else 'utf-8-sig')
    assert entry['raw_bytes'] == len(raw)
    assert entry['normalized_bytes'] == len(text.encode('utf-8'))
    assert entry['normalized_sha256'] == hashlib.sha256(text.encode('utf-8')).hexdigest()
    write_source(tmp_path, raw=text.encode(encoding))
    assert source_identity(require_sources(tmp_path, manifest)) == source_identity(manifest)


def test_empty_python_source_is_valid(tmp_path):
    write_source(tmp_path, raw=b'')
    manifest = capture_sources(tmp_path, [PATH])
    entry = manifest['files'][0]
    assert entry['raw_bytes'] == entry['normalized_bytes'] == 0
    assert entry['raw_sha256'] == entry['normalized_sha256'] == hashlib.sha256(b'').hexdigest()
    assert require_sources(tmp_path, manifest) == manifest


def test_noncanonical_encoding_cannot_hide_source_edits(tmp_path):
    # Python accepts UTF-7's alternate spelling of ASCII, which decodes away.
    raw = b'# coding: utf-7\nVALUE = +ADE-\n'
    compile(raw, '<fixture>', 'exec')
    assert raw.decode('utf-7').endswith('VALUE = 1\n')
    write_source(tmp_path, raw=raw)
    with pytest.raises(ValueError, match='invalid Python source or encoding'):
        capture_sources(tmp_path, [PATH])


def test_capture_compiles_without_executing_or_importing_sources(tmp_path):
    write_source(tmp_path, raw=b'import nonexistent_module\nraise RuntimeError("never run")\n')
    manifest = capture_sources(tmp_path, [PATH])
    assert require_sources(tmp_path, manifest) == manifest


def test_module_loads_independently_of_gabes_imports(monkeypatch):
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == 'gabes' or name.startswith('gabes.'):
            raise AssertionError('provenance must not import the mutable gabes graph')
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', guarded_import)
    namespace = runpy.run_path(provenance.__file__)
    assert all(callable(namespace[name]) for name in provenance.__all__)


@pytest.mark.parametrize('field,value', [
    ('schema', 'legacy'), ('schema', 1), ('source_sha256', '0'*64),
    ('manifest_sha256', '0'*64), ('source_sha256', 'X'*64),
    ('manifest_sha256', None), ('files', []), ('files', {}), ('files', [None]),
])
def test_manifest_schema_shape_and_digest_corruption_is_rejected(captured, field, value):
    root, _, manifest = captured
    manifest[field] = value
    with pytest.raises(ValueError):
        source_identity(manifest)
    with pytest.raises(ValueError):
        require_sources(root, manifest)


@pytest.mark.parametrize('field,value', [
    ('path', '../escape.py'), ('path', 'renamed.py'), ('path', 'model.txt'),
    ('raw_sha256', '0'*64), ('normalized_sha256', '0'*64),
    ('raw_sha256', 'G'*64), ('normalized_sha256', 'abc'),
    ('raw_bytes', 999), ('raw_bytes', -1), ('raw_bytes', True),
    ('normalized_bytes', 999), ('normalized_bytes', 1.0),
    ('encoding', 'latin-1'), ('encoding', 'unknown-codec'),
    ('encoding', ''), ('encoding', None), ('encoding', 'base64'),
])
def test_record_metadata_and_digest_tampering_is_rejected(captured, field, value):
    root, _, manifest = captured
    manifest['files'][0][field] = value
    with pytest.raises(ValueError):
        source_identity(manifest)
    with pytest.raises(ValueError):
        require_sources(root, manifest)


def test_manifest_missing_extra_duplicate_and_unsorted_records_are_rejected(captured):
    root, _, manifest = captured
    write_source(root, 'second.py')
    both = capture_sources(root, [PATH, 'second.py'])
    variants = [None, [], {PATH: manifest['files'][0]['raw_sha256']}]
    for key in manifest:
        bad = deepcopy(manifest)
        del bad[key]
        variants.append(bad)
    for key in manifest['files'][0]:
        bad = deepcopy(manifest)
        del bad['files'][0][key]
        variants.append(bad)
    variants.append(dict(manifest, extra=True))
    bad = deepcopy(manifest)
    bad['files'][0]['extra'] = True
    variants.append(bad)
    bad = deepcopy(manifest)
    bad['files'].append(deepcopy(bad['files'][0]))
    variants.append(bad)
    bad = deepcopy(both)
    bad['files'].reverse()
    variants.append(bad)
    bad = deepcopy(both)
    bad['files'].pop()
    variants.append(bad)
    for bad in variants:
        with pytest.raises(ValueError):
            source_identity(bad)
        with pytest.raises(ValueError):
            require_sources(root, bad)
