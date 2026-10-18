"""Semantic corruption controls, confined to fresh temporary copies."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
FLAGS = ['-B'] + (['-' + 'O' * sys.flags.optimize] if sys.flags.optimize else [])

def trial(name, mutation, command, expected):
    with tempfile.TemporaryDirectory(prefix='dormancy-negative-') as folder:
        target = Path(folder)
        for path in ROOT.iterdir():
            if path.suffix in {'.py', '.js', '.json'}:
                shutil.copy2(path, target / path.name)
        mutation(target)
        result = subprocess.run([sys.executable, *FLAGS, *command], cwd=target,
                                capture_output=True, text=True, timeout=180)
        if result.returncode == 0 or expected not in result.stderr:
            raise RuntimeError(f'{name}: expected semantic rejection absent: {result.stderr}')
        return {'control': name, 'rejected': True, 'exit_code': result.returncode,
                'failure': result.stderr.splitlines()[-1]}

def polynomial(target):
    path = target / 'hopf_polynomial_coefficients.json'
    values = json.loads(path.read_text())
    values[0] = str(int(values[0]) + 1)
    path.write_text(json.dumps(values))

def vector_field(target):
    path = target / 'model.py'
    source = path.read_text()
    old = "DV*a-(r+v)*i,"
    if source.count(old) != 1:
        raise RuntimeError('mutation target not unique')
    path.write_text(source.replace(old, 'DV*a-(r+v)*i+1,'))

if __name__ == '__main__':
    results = [
        trial('corrupted Hopf constant coefficient', polynomial, ['verify_symbolic.py'], 'AssertionError'),
        trial('altered infection balance', vector_field, ['verify_symbolic.py'], 'AssertionError'),
        trial('inverted interval endpoints', lambda target: None,
              ['-c', 'from interval_arithmetic import I; I.raw(2,1)'], 'AssertionError')]
    print(json.dumps({'status': 'PASS', 'optimized': sys.flags.optimize, 'controls': results}, indent=2))
