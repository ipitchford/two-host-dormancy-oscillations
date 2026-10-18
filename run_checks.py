from check_support import require
"""Portable, fail-fast reproduction of the dormancy proof certificates.

Run from any working directory with Python 3.10+, Node 18+, and the pinned
SymPy dependencies installed. Outputs are local to this paper directory.
"""
import datetime,hashlib,importlib.metadata,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/'reproduction'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    node=shutil.which('node') or os.environ.get('CODEX_PRIMARY_RUNTIME_NODE')
    if not node:raise SystemExit('Node.js 18 or later is required for the four recovered exact checks.')
    for pkg,want in [('sympy','1.14.0'),('mpmath','1.3.0')]:
        actual=importlib.metadata.version(pkg)
        if actual!=want:raise SystemExit(f'{pkg} {want} required; found {actual}')
    scripts=[('exact_counterexample.js',node),('all_c_hurwitz.js',node),
             ('hopf_polynomial.js',node),('q_coefficients_samples.js',node),
             ('interval_arithmetic.py',sys.executable),
             ('verify_symbolic.py',sys.executable),('certify_hopf.py',sys.executable)]
    inputs=[p for p in HERE.iterdir() if p.suffix in ['.py','.js']]+[HERE/'hopf_polynomial_coefficients.json']
    before={p.name:sha(p) for p in sorted(inputs)}
    record={'status':'RUNNING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'working_directory':str(HERE),'versions':{'python':platform.python_version(),
        'node':subprocess.check_output([node,'--version'],text=True).strip(),
        'sympy':importlib.metadata.version('sympy'),'mpmath':importlib.metadata.version('mpmath')},
      'environment':{'PYTHONPATH':os.environ.get('PYTHONPATH','')},'input_sha256':before,'checks':[]}
    start=time.monotonic()
    for name,exe in scripts:
        print('Running '+name,flush=True);t=time.monotonic()
        flags = ['-B'] + (['-' + 'O' * sys.flags.optimize] if sys.flags.optimize else [])
        command=[exe]+(flags if name.endswith('.py') else [])+[str(HERE/name)]
        done=subprocess.run(command,cwd=HERE,text=True,capture_output=True,timeout=300)
        elapsed=time.monotonic()-t;stem=Path(name).stem
        (OUT/(stem+'.stdout.txt')).write_text(done.stdout)
        (OUT/(stem+'.stderr.txt')).write_text(done.stderr)
        item={'script':name,'command':command,'exit_code':done.returncode,'elapsed_seconds':elapsed,
          'status':'passed' if done.returncode==0 else 'failed',
          'stdout':stem+'.stdout.txt','stderr':stem+'.stderr.txt'}
        record['checks'].append(item)
        if done.returncode:
            record['status']='FAIL';record['elapsed_seconds']=time.monotonic()-start
            (OUT/'run_record.json').write_text(json.dumps(record,indent=2)+'\n')
            print(done.stderr,file=sys.stderr);raise SystemExit(done.returncode)
        if name.endswith('.js'):
            result=json.loads(done.stdout);require(result['status']=='PASS', "check failed: result['status']=='PASS'")
            (OUT/(stem+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    # Concrete counts and equation fingerprints, not just program exit codes.
    j=lambda n:json.loads((OUT/(n+'.json')).read_text())
    require(len(j('exact_counterexample')['exactInterpolationPoints'])==7, "check failed: len(j('exact_counterexample')['exactInterpolationPoints'])==7")
    require(j('exact_counterexample')['coexistenceRouthSigns']==list('+++++-+'), "check failed: j('exact_counterexample')['coexistenceRouthSigns']==list('+++++-+')")
    require(j('exact_counterexample')['trace']=='-2291/100', "check failed: j('exact_counterexample')['trace']=='-2291/100'")
    require([x['sign'] for x in j('all_c_hurwitz')['hurwitz']]==list('++++--'), "check failed: [x['sign'] for x in j('all_c_hurwitz')['hurwitz']]==list('++++--')")
    require(len(j('q_coefficients_samples')['q_values'])==13, "check failed: len(j('q_coefficients_samples')['q_values'])==13")
    require(j('hopf_polynomial')['coefficients_in_ascending_order']==json.loads((HERE/'hopf_polynomial_coefficients.json').read_text()), "check failed: j('hopf_polynomial')['coefficients_in_ascending_order']==json.loads((HERE/'hopf_polynomial_coefficients.json').read_text())")
    hopf=json.loads((HERE/'hopf_interval_certificate.json').read_text())
    require(hopf['status']=='PASS' and hopf['decimal_grid_digits']==60, "check failed: hopf['status']=='PASS' and hopf['decimal_grid_digits']==60")
    require(before=={p.name:sha(p) for p in sorted(inputs)}, 'checker sources changed during execution')
    record['status']='PASS';record['elapsed_seconds']=time.monotonic()-start
    record['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    record['counts']={'programs_passed':7,'exact_C_interpolation_points':7,
      'all_C_Hurwitz_signs':'++++--','q_sample_points':13,'Hopf_polynomial_degree':17,
      'exact_rational_arithmetic_tests':99,'interval_decimal_grid_digits':60}
    record['output_sha256']={n:sha(HERE/n) for n in ['family_symbolic_certificate.json','hopf_symbolic_certificate.json','hopf_interval_certificate.json']}
    (OUT/'run_record.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'status':record['status'],'elapsed_seconds':record['elapsed_seconds'],'counts':record['counts']},indent=2))


if __name__=='__main__':main()
