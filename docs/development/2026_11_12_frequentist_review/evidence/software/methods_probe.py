"""Fresh-process checks of methods provenance, only temporary workspace mutations."""
import argparse,csv,hashlib,json,pathlib,shutil,subprocess,sys,tempfile
ap=argparse.ArgumentParser();ap.add_argument('repo',type=pathlib.Path);a=ap.parse_args();source=a.repo.resolve()
with tempfile.TemporaryDirectory(dir='.') as tmp:
    root=pathlib.Path(tmp).resolve()/'app';root.mkdir()
    for d in ['data','tools/early_tox','fit/stan','draws','R']:
        shutil.copytree(source/d,root/d,ignore=shutil.ignore_patterns('__pycache__','*.csv.gz'))
    unit='docs/development/2026_11_10_early_toxicity/step5_review/overdispersed/U02';fit='fit/output_refit_2026_10_25c_main'
    for f in ['fit/fit_refit.py','fit/gate.py','tools/build_model_methods.py',f'{fit}/di_primary/identity.json',f'{fit}/dd_calibrated_onset7_n1000/identity.json',f'{unit}/manifest.json',f'{unit}/diagnostics.json']:
        dest=root/f;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source/f,dest)
    def run(label):
        p=subprocess.run([sys.executable,'-B',str(root/'tools/build_model_methods.py'),'--root',str(root),'--check'],capture_output=True,text=True)
        print(label,'exit',p.returncode,(p.stdout+p.stderr).strip()[:300]);return p.returncode
    run('unchanged positive')
    for name,column in [('tumor_di_primary.csv','alpha'),('early_fit_draws.csv','r_e')]:
        path=root/'draws'/name;original=path.read_bytes();lines=original.decode().splitlines();header=lines[0].split(',');values=lines[1].split(',');values[header.index(column)]='100';lines[1]=','.join(values);path.write_text('\n'.join(lines)+'\n')
        print('changed',name,'column',column,'file hash altered',hashlib.sha256(path.read_bytes()).hexdigest()!=hashlib.sha256(original).hexdigest())
        run('changed draw bytes methods check');path.write_bytes(original)
    stan=root/'tools/early_tox/step5_model.stan';original=stan.read_text();stan.write_text(original.replace('sigma_arm ~ normal(0, od_sd);','sigma_arm ~ normal(0, 2 * od_sd);'))
    run('changed scientific model negative')
