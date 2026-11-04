"""Fetch pinned third-party assets for local and hosted documentation builds."""
import hashlib
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen


CACHE = Path(__file__).resolve().parent / '_build' / 'assets'
MATHJAX = (
    ('mathjax-3.2.2.js',
     'https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg-full.js',
     'a4354ff94fd868aea0cc6eaaa79a57fda0588646fc46ee3700a349ee0a11cbe6'),
    ('mathjax-LICENSE',
     'https://cdn.jsdelivr.net/npm/mathjax@3.2.2/LICENSE',
     'cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30'),
)


def install_mathjax(destination):
    """Verify cached or downloaded assets, then copy them into the site."""
    CACHE.mkdir(parents=True, exist_ok=True)
    destination.mkdir(parents=True, exist_ok=True)
    for name, url, checksum in MATHJAX:
        cached = CACHE / name
        if cached.is_file():
            data = cached.read_bytes()
        else:
            with urlopen(url, timeout=60) as response:
                data = response.read()
        if hashlib.sha256(data).hexdigest() != checksum:
            raise ValueError(f'MathJax checksum mismatch for {name}; '
                             f'remove {cached} if cached and retry the build.')
        if not cached.is_file():
            # Atomic replacement also permits concurrent documentation builds.
            with tempfile.NamedTemporaryFile(dir=CACHE, delete=False) as stream:
                temporary = Path(stream.name)
                try:
                    stream.write(data)
                    stream.close()
                    temporary.replace(cached)
                finally:
                    temporary.unlink(missing_ok=True)
        shutil.copyfile(cached, destination / name)
