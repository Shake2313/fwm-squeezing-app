"""Portable source manifests for new runs; legacy hashes are never upgraded.

``capture_sources(root, paths)`` returns a JSON-compatible manifest with sorted
``files`` records: path, detected encoding, raw/normalized SHA-256 and byte counts.
Normalized bytes are decoded Python source encoded as UTF-8 after CRLF/CR -> LF.
``source_sha256`` binds the schema, filenames, encoding and normalized contents;
``manifest_sha256`` additionally binds the raw audit metadata. Neither contains
the checkout root. ``source_identity`` validates both checksums before returning
the former; ``require_sources`` checks the files and returns a fresh manifest.

Checksums detect corruption, not authenticity: someone who replaces a manifest
and recomputes its checksums can forge it. Callers must retain a trusted manifest
and explicitly enumerate every source they depend on; imports are not followed.
Disk hashes do not establish the origin of code already loaded in an interpreter.
Snapshots do not lock files against concurrent edits during or after capture.
"""

import codecs
import hashlib
import json
from pathlib import Path
import re
import tokenize


__all__ = ['capture_sources', 'source_identity', 'require_sources']

SCHEMA = 'gabes.python-source-manifest.v1'
_MANIFEST_KEYS = {'schema', 'files', 'source_sha256', 'manifest_sha256'}
_FILE_KEYS = {'path', 'encoding', 'raw_sha256', 'raw_bytes',
              'normalized_sha256', 'normalized_bytes'}
_IDENTITY_KEYS = ('path', 'encoding', 'normalized_sha256', 'normalized_bytes')
_SHA256 = re.compile(r'[0-9a-f]{64}')
_WINDOWS_DEVICE = re.compile(r'(CON|PRN|AUX|NUL|COM[1-9¹²³]|LPT[1-9¹²³])', re.I)


def _digest(value):
    data = json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=True, allow_nan=False).encode('ascii')
    return hashlib.sha256(data).hexdigest()


def _path_name(value):
    if type(value) is not str or not value or not value.endswith('.py'):
        raise ValueError('source paths must be explicit repo-relative POSIX .py strings')
    parts = value.split('/')
    if (any(part in ('', '.', '..') or part.endswith((' ', '.')) for part in parts)
            or any(char in value for char in '<>:"\\|?*')
            or any(ord(char) < 32 or 0xD800 <= ord(char) <= 0xDFFF for char in value)
            or any(_WINDOWS_DEVICE.fullmatch(part.split('.')[0]) for part in parts)):
        raise ValueError(f'unsafe or non-portable source path: {value!r}')
    return value


def _paths(paths):
    if isinstance(paths, (str, bytes, dict)):
        raise ValueError('paths must be a nonempty iterable of explicit source paths')
    try:
        names = [_path_name(path) for path in paths]
    except TypeError as exc:
        raise ValueError('paths must be a nonempty iterable of explicit source paths') from exc
    if not names:
        raise ValueError('at least one source path is required')
    # Case aliases are ambiguous when a checkout moves to Windows or macOS.
    if len({name.casefold() for name in names}) != len(names):
        raise ValueError('duplicate source paths (including case aliases)')
    return sorted(names)


def _source_digest(manifest):
    return _digest({'schema': manifest['schema'], 'files': [
        {key: entry[key] for key in _IDENTITY_KEYS} for entry in manifest['files']]})


def _manifest_digest(manifest):
    return _digest({key: value for key, value in manifest.items()
                    if key != 'manifest_sha256'})


def _check_sha256(value):
    if type(value) is not str or _SHA256.fullmatch(value) is None:
        raise ValueError('SHA-256 digests must be 64 lowercase hexadecimal characters')


def source_identity(manifest) -> str:
    """Validate a v1 manifest and return its portable source SHA-256 identity.

    This checks internal integrity without reading files. Missing, extra or
    legacy fields are rejected, including edits to the raw audit metadata.
    """
    if type(manifest) is not dict or set(manifest) != _MANIFEST_KEYS:
        raise ValueError('invalid source manifest shape')
    if type(manifest['schema']) is not str or manifest['schema'] != SCHEMA:
        raise ValueError('unsupported source manifest schema; capture new runs only')
    files = manifest['files']
    if type(files) is not list or not files:
        raise ValueError('manifest files must be a nonempty list')
    for entry in files:
        if type(entry) is not dict or set(entry) != _FILE_KEYS:
            raise ValueError('invalid source file record shape')
        _path_name(entry['path'])
        for key in ('raw_sha256', 'normalized_sha256'):
            _check_sha256(entry[key])
        for key in ('raw_bytes', 'normalized_bytes'):
            if type(entry[key]) is not int or entry[key] < 0:
                raise ValueError(f'{key} must be a nonnegative integer')
        encoding = entry['encoding']
        if type(encoding) is not str or not encoding:
            raise ValueError('source encoding must be a known text encoding')
        try:
            codecs.lookup(encoding)
            b''.decode(encoding)
        except (LookupError, UnicodeError, ValueError) as exc:
            raise ValueError('source encoding must be a known text encoding') from exc
    names = [entry['path'] for entry in files]
    if names != _paths(names):
        raise ValueError('manifest files must be sorted by relative path')
    _check_sha256(manifest['source_sha256'])
    _check_sha256(manifest['manifest_sha256'])
    if manifest['source_sha256'] != _source_digest(manifest):
        raise ValueError('source manifest identity digest mismatch')
    if manifest['manifest_sha256'] != _manifest_digest(manifest):
        raise ValueError('source manifest audit digest mismatch')
    return manifest['source_sha256']


def capture_sources(root, paths) -> dict:
    """Capture explicit Python sources below root without importing/executing them.

    Paths must be portable relative POSIX strings, with no duplicates or aliases
    to the same file. Symlinks within root are allowed; escapes are rejected.
    Source encoding must round-trip losslessly without alternative byte spellings.
    All validation, missing-file and invalid-Python failures raise ValueError.
    """
    names = _paths(paths)
    try:
        base = Path(root).resolve(strict=True)
        if not base.is_dir():
            raise ValueError('source root must be a directory')
    except (OSError, RuntimeError, TypeError) as exc:
        raise ValueError('source root must be an existing directory') from exc
    files = []
    seen = set()
    for name in names:
        try:
            path = base.joinpath(*name.split('/')).resolve(strict=True)
            if not path.is_relative_to(base):
                raise ValueError(f'source path escapes root: {name}')
            if not path.is_file():
                raise ValueError(f'source is not a regular file: {name}')
            stat = path.stat()
            file_id = (stat.st_dev, stat.st_ino)
            if file_id in seen:
                raise ValueError(f'duplicate source file alias: {name}')
            seen.add(file_id)
            raw = path.read_bytes()
        except (OSError, RuntimeError) as exc:
            raise ValueError(f'cannot read source file: {name}') from exc
        try:
            # bytes.splitlines recognizes only CR, LF and CRLF. Detection must
            # inspect physical lines even in a CR-only Python source file.
            lines = iter(raw.splitlines(keepends=True))
            encoding, _ = tokenize.detect_encoding(lambda: next(lines, b''))
            decoded = raw.decode(encoding)
            # A codec must not collapse distinct byte spellings of source text.
            if decoded.encode(encoding) != raw:
                raise ValueError('source encoding does not round-trip losslessly')
            # EOL equivalence is Python universal newline semantics only.
            # There is no reduction of numerical/physical validation.
            # Preserve all strings, comments and whitespace; never compare ASTs.
            normalized = decoded.replace('\r\n', '\n').replace('\r', '\n')
            compile(raw, name, 'exec', dont_inherit=True, optimize=0)
            compile(normalized, name, 'exec', dont_inherit=True, optimize=0)
            normalized_bytes = normalized.encode('utf-8')
        except (SyntaxError, UnicodeError, LookupError, ValueError) as exc:
            raise ValueError(f'invalid Python source or encoding: {name}') from exc
        files.append({'path': name, 'encoding': encoding,
                      'raw_sha256': hashlib.sha256(raw).hexdigest(),
                      'raw_bytes': len(raw),
                      'normalized_sha256': hashlib.sha256(normalized_bytes).hexdigest(),
                      'normalized_bytes': len(normalized_bytes)})
    manifest = {'schema': SCHEMA, 'files': files}
    manifest['source_sha256'] = _source_digest(manifest)
    manifest['manifest_sha256'] = _manifest_digest(manifest)
    return manifest


def require_sources(root, manifest) -> dict:
    """Reject source changes and return the current raw snapshot for audit.

    Only the same explicit paths are recaptured. Relocation and Python universal
    newlines may change raw audit fields, but not the verified source identity.
    The supplied manifest is never mutated or relabelled.
    """
    expected = source_identity(manifest)
    current = capture_sources(root, [entry['path'] for entry in manifest['files']])
    if current['source_sha256'] != expected:
        changed = [old['path'] for old, new in zip(manifest['files'], current['files'])
                   if any(old[key] != new[key] for key in _IDENTITY_KEYS)]
        raise ValueError('Python source changed: ' + ', '.join(changed))
    return current
