"""Regression cases: reject invalid values and preserve the minimal input contract."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from prepare_gse_inputs import write_xlsx, prepare
from validate_promotheus_inputs import first_sheet_rows, validate_workbook

class StrictValidationTests(unittest.TestCase):
    def test_minimal_input_gets_numeric_gse(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)
            for kind, accession in [('Microarray','GSE000001'),('RNA_Seq','GSE000002')]:
                f=p/(kind+'.xlsx')
                write_xlsx(f,[['Accession','taxon','reference_genome'],[accession,'Toy','ASM100001v1']])
                self.assertEqual(validate_workbook(f)[1],[])
            prepare(p,p/'shards')
            for f in (p/'shards').glob('*/*.xlsx'):
                rows=first_sheet_rows(f)
                self.assertEqual(rows[1][rows[0].index('GSE')],rows[1][rows[0].index('Accession')][3:])

    def test_output_content(self):
        original=first_sheet_rows(ROOT/'tests/fixtures/toy_expected/final_table.xlsx')
        for case in ('valid','nan','infinity','empty_numeric','duplicate_gene','empty_gene','tsv_value','duplicate_column','missing_percentile'):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temp:
                p=Path(temp);(p/'results').mkdir()
                rows=[r[:] for r in original]
                if case in ('nan','infinity','empty_numeric'):
                    rows[1][rows[0].index('Mean_Expression')]={'nan':'NaN','infinity':'Inf','empty_numeric':''}[case]
                if case=='duplicate_gene':rows[2][rows[0].index('GeneID')]=rows[1][rows[0].index('GeneID')]
                if case=='empty_gene':rows[1][rows[0].index('GeneID')]=''
                if case=='duplicate_column':
                    for r in rows:r.append(r[2])
                if case=='missing_percentile':
                    idx=rows[0].index('Percentile_GSE000001')
                    for r in rows:r.pop(idx)
                write_xlsx(p/'results/final_table.xlsx',rows)
                tsv=[r[:] for r in rows]
                if case=='tsv_value':tsv[1][tsv[0].index('Mean_Expression')]='999999'
                with (p/'results/final_table.tsv').open('w') as f:csv.writer(f,delimiter='\t').writerows(tsv)
                (p/'orthologous_genes.txt').write_text('toy\n')
                (p/'record.tsv').write_text('field\tvalue\nexport_tsv\ttrue\n')
                (p/'input.tsv').write_text('accessions\nGSE000001,GSE000002\n')
                (p/'check.tsv').write_text('accession\tstatus\nGSE000001\tcomplete\nGSE000002\tcomplete\n')
                cmd=[sys.executable,str(ROOT/'bin/verify_promotheus_outputs.py'),'--output-dir',str(p),'--run-record',str(p/'record.tsv'),'--input-manifest',str(p/'input.tsv'),'--checkpoint-manifest',str(p/'check.tsv'),'--manifest',str(p/'manifest.tsv'),'--report',str(p/'report.json')]
                result=subprocess.run(cmd,capture_output=True,text=True)
                self.assertEqual(result.returncode,0 if case=='valid' else 1,result.stdout+result.stderr)
                self.assertEqual(json.loads((p/'report.json').read_text())['status'],'pass' if case=='valid' else 'fail')

if __name__=='__main__':unittest.main()
