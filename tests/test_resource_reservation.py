#!/usr/bin/env python3
"""Check balanced scheduling with real Nextflow and deterministic toy work."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ResourceReservationTests(unittest.TestCase):
    def test_reservation_overlap_resume_and_insufficient_budget(self):
        with tempfile.TemporaryDirectory(prefix='promotheus-balanced-') as temp:
            root = Path(temp)
            project = root / 'project'
            project.mkdir()
            for name in ('main.nf', 'nextflow.config'):
                shutil.copy2(ROOT / name, project)
            for name in ('bin', 'Scripts', 'tests/fixtures'):
                shutil.copytree(ROOT / name, project / name, ignore=shutil.ignore_patterns('__pycache__'))
            runner = project / 'bin/run_modular_stage.py'
            source = runner.read_text()
            self.assertEqual(source.count('range(1, 4)'), 1)
            source = source.replace('range(1, 4)', 'range(1, 13)')
            source = source.replace('time.sleep(float(config.get("mock_delay", 0)))',
                                    'time.sleep({"orthology_pair": 3, "rna_download": 1}.get(stage, 0.1))')
            runner.write_text(source)
            base = ['nextflow', 'run', str(project / 'main.nf'), '-profile', 'test,balanced',
                    '--input_dir', str(project / 'tests/fixtures/toy_input'),
                    '--output_dir', str(root / 'output'), '--asmcode', 'ASM317683v1',
                    '--threads', '8', '--memory', '4096', '--task_cpus', '1', '--task_memory', '2048',
                    '--orthology_cpus', '1', '--orthology_memory', '512',
                    '--max_orthology_parallel', '6', '-ansi-log', 'false']
            def execute(name, options=(), success=True):
                trace = root / (name + '.tsv')
                result = subprocess.run(base + ['-with-trace', str(trace)] + list(options), cwd=root,
                                        env=dict(os.environ, NXF_OFFLINE='true'), text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
                self.assertEqual(result.returncode == 0, success, result.stdout)
                rows = list(csv.DictReader(trace.read_text().splitlines(), delimiter='\t')) if trace.exists() else []
                return result.stdout, rows
            log, rows = execute('balanced')
            self.assertIn('concurrency capped at 4 (requested 6)', log)
            records = []
            for row in rows:
                work = Path(row['workdir'])
                if (work / 'task.json').exists():
                    config = json.loads((work / 'task.json').read_text())
                    records.append(json.loads((work / config['output'] / 'task_record.json').read_text()))
            pairs = [r for r in records if r['stage'] == 'orthology_pair']
            self.assertEqual(len(pairs), 12)
            peak = max(sum(p['started'] <= t['started'] < p['finished'] for p in pairs) for t in pairs)
            self.assertGreaterEqual(peak, 2)
            self.assertLessEqual(peak, 4)
            align = next(r for r in records if r['stage'] == 'rna_align')
            self.assertLess(align['started'], max(p['finished'] for p in pairs))
            download = next(r for r in records if r['stage'] == 'rna_download')
            print(f"PASS: peak orthology concurrency {peak}; RNA alignment wait {align['started'] - download['finished']:.3f}s; alignment overlaps orthology")
            _, cached = execute('resume', ['-resume', '--reserve_rnaseq_resources', 'false'])
            self.assertTrue(all(r['status'] == 'CACHED' for r in cached), cached)
            log, _ = execute('insufficient', ['--memory', '2048'], success=False)
            self.assertIn('Cannot reserve one RNA-Seq task', log)


if __name__ == '__main__':
    unittest.main()
