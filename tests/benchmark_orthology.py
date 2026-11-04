#!/usr/bin/env python3
"""Compare fresh serial/parallel orthology runs using the same local genomes."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--genomes', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--asmcode', required=True)
    args = parser.parse_args()
    project, inputs, genomes, root = [p.resolve() for p in (args.project, args.input, args.genomes, args.output)]
    root.mkdir(parents=True, exist_ok=False)
    def state(status, **details):
        (root / 'status.json').write_text(json.dumps(dict(status=status, updated_at=datetime.now(timezone.utc).isoformat(), **details), indent=2)+'\n')
    def snapshot():
        files = [project/'main.nf', project/'nextflow.config']
        for folder in (project/'bin', project/'Scripts', project/'tests/fixtures', inputs, genomes):
            files.extend(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc')
        return {str(p):sha(p) for p in sorted(set(files))}
    frozen = snapshot()
    (root/'snapshot_sha256.json').write_text(json.dumps(frozen,indent=2)+'\n')
    reports = {}
    try:
        for mode, forks in [('serial',1), ('parallel',2)]:
            state('running', mode=mode)
            work = root/mode
            work.mkdir()
            cmd = ['nextflow','run',str(project/'main.nf'),'-entry','ORTHOLOGY_ONLY','-profile','standard',
                   '--input_dir',str(inputs),'--output_dir',str(work/'output'),'--orthology_genomes',str(genomes),
                   '--asmcode',args.asmcode,'--threads','4','--memory','12000','--task_cpus','1','--task_memory','3000',
                   '--orthology_cpus','1','--orthology_memory','3000','--max_orthology_parallel',str(forks),
                   '-work-dir',str(work/'work'),'-ansi-log','false']
            (work/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
            started = time.monotonic()
            with (work/'execution.log').open('w') as log:
                subprocess.run(cmd,cwd=work,stdout=log,stderr=subprocess.STDOUT,check=True,
                               env=dict(os.environ,NXF_OFFLINE='true'))
            elapsed = time.monotonic()-started
            traces = list((work/'output/pipeline_info').glob('execution_trace_*.txt'))
            assert len(traces)==1, 'Expected one trace per fresh run'
            tasks = list(csv.DictReader(traces[0].read_text().splitlines(),delimiter='\t'))
            assert tasks and all(t['status']=='COMPLETED' and t['exit']=='0' for t in tasks), 'Failed/cached tasks'
            pairs = [t for t in tasks if 'ORTHOLOGY_PAIR' in t['name']]
            assert len(pairs)>=2 and len(tasks)==len(pairs)+3, 'Unexpected task coverage'
            stamp = datetime.fromisoformat
            events=sorted((stamp(t[k]),delta) for t in pairs for k,delta in [('start',1),('complete',-1)])
            running=peak=0
            for _,delta in events:
                running+=delta; peak=max(peak,running)
            assert peak==forks, 'Requested concurrency was not achieved'
            out=work/'output/orthology/orthologous_genes.txt'
            mapping_hashes={}
            for task in pairs:
                config=json.loads((Path(task['workdir'])/'task.json').read_text())
                mapping=Path(task['workdir'])/config['output']/'mapping.rds'
                mapping_hashes[config['key']]=sha(mapping)
            reports[mode]={'command_wall_seconds':elapsed,'tasks':len(tasks),'pairs':len(pairs),
                           'workflow_seconds':(max(stamp(t['complete']) for t in tasks)-min(stamp(t['start']) for t in tasks)).total_seconds(),
                           'pair_window_seconds':(max(stamp(t['complete']) for t in pairs)-min(stamp(t['start']) for t in pairs)).total_seconds(),
                           'peak_concurrent_pairs':peak,'output_sha256':sha(out),'mapping_sha256':mapping_hashes,
                           'trace':str(traces[0])}
            (root/'progress.json').write_text(json.dumps(reports,indent=2)+'\n')
        assert snapshot()==frozen, 'Sources or inputs changed during the comparison'
        assert reports['serial']['output_sha256']==reports['parallel']['output_sha256'], 'Final orthology differs'
        assert reports['serial']['mapping_sha256']==reports['parallel']['mapping_sha256'], 'Pair mappings differ'
        report={'status':'pass','runs':reports,'identical_outputs':True,
                'workflow_speedup':reports['serial']['workflow_seconds']/reports['parallel']['workflow_seconds'],
                'limits':'One serial/parallel pair on local archived genomes; no Nextflow resume. OS/Docker filesystem caches are not flushed. Fixed serial-first order and host load can affect timing. Orthology-only measurement, not end-to-end speedup; CPU and RSS from nested Docker are not inferred from the Nextflow trace.'}
        (root/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
        state('complete')
        print(json.dumps(report,indent=2))
    except BaseException as error:
        state('failed',error=str(error))
        raise


if __name__=='__main__':
    main()
