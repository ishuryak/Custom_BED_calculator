import argparse, hashlib, json, pathlib, subprocess
p=argparse.ArgumentParser();p.add_argument('original',type=pathlib.Path);p.add_argument('snapshot',type=pathlib.Path);p.add_argument('gitdir');p.add_argument('revision');a=p.parse_args()
git=['C:/Program Files/Git/cmd/git.exe','--git-dir='+a.gitdir]
tree=subprocess.check_output(git+['ls-tree','-r',a.revision],text=True).splitlines()
records=[]
for line in tree:
    info,path=line.split('\t',1);expected=info.split()[2]
    raw=(a.snapshot/path).read_bytes(); original=(a.original/path).read_bytes() if (a.original/path).exists() else None
    blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    records.append(dict(path=path,git_blob=expected,snapshot_blob=blob,sha256=hashlib.sha256(raw).hexdigest(),original_equal=raw==original))
pathlib.Path('source_manifest.json').write_text(json.dumps(records,indent=1))
print('tracked',len(records),'snapshot blob mismatches',[r['path'] for r in records if r['git_blob']!=r['snapshot_blob']])
print('working-tree differences',[r['path'] for r in records if not r['original_equal']])
