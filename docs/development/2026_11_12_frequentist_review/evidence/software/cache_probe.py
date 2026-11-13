"""Current cache byte-binding regression with one copied real cache; no fitting."""
import argparse,importlib.util,json,pathlib,shutil,sys,tempfile,types
import numpy as np
p=argparse.ArgumentParser();p.add_argument('repo',type=pathlib.Path);p.add_argument('cache',type=pathlib.Path);a=p.parse_args();root=a.repo.resolve()
def deny(*a,**k):raise RuntimeError('review never samples')
st=types.ModuleType('cmdstanpy');st.__version__='1.3.0';st.CmdStanModel=deny;st.cmdstan_path=lambda:'/stub/cmdstan-2.40.0';sys.modules['cmdstanpy']=st
spec=importlib.util.spec_from_file_location('s6',root/'tools/early_tox/step6_validate.py');s6=importlib.util.module_from_spec(spec);spec.loader.exec_module(s6)
rows,fams=s6.s5.load_rows();name='candidate__all_primary';settings=s6.MODELS['candidate']
version=np.__version__;np.__version__='2.4.3' # Recorded fit environment in build_model_methods.EARLY_ENV, for identity reconstruction only.
try:
    with tempfile.TemporaryDirectory(dir='.') as td:
        target=pathlib.Path(td);file=target/(name+'.npz');shutil.copyfile(a.cache/file.name,file);shutil.copyfile(a.cache/'cache_manifest.json',target/'cache_manifest.json')
        x=s6.fit_cached(None,name,rows,fams,settings,str(target),rescore_only=True)
        print('positive real cache identity reproduced and loaded',len(x['r']),'draws')
        x['r'][0]=150.;np.savez_compressed(file,**x)
        try:s6.fit_cached(None,name,rows,fams,settings,str(target),rescore_only=True)
        except RuntimeError as e:print('old changed-array defect refused',str(e).replace(str(target),'TEMP'))
finally:np.__version__=version
print('Native NumPy',version,'used for arrays; fit metadata version used only for deterministic hash reconstruction')
