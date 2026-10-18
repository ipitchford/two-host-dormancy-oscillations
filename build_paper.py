"""Rebuild the paper from shipped Markdown using Pandoc and Tectonic."""
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
title, body = (ROOT / 'manuscript.md').read_text().split('\n', 1)
output = ROOT / 'paper-build'
output.mkdir(exist_ok=True)
source = output / 'body.md'
# The generated TeX lives one directory below the shipped figures.
source.write_text(body)
tex = output / 'paper.tex'
subprocess.run(['pandoc', str(source), '--standalone', '--shift-heading-level-by=-1',
                '--from=markdown+tex_math_dollars+autolink_bare_uris', '--to=latex',
                '-M', 'title='+title.removeprefix('# '), '-M', 'author=Anonymous',
                '-M', 'date=18 October 2026', '-V', 'geometry:margin=25mm',
                '-V', 'fontsize=11pt', '-V', 'colorlinks=true',
                '-V', 'documentclass=article', '-o', str(tex)], check=True)
subprocess.run(['tectonic', '--untrusted', '--keep-logs', 'paper.tex'],
               cwd=output, check=True)
