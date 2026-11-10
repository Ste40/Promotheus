from pathlib import Path
from html import escape
W,H=2000,820
s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="title desc">', '<title id="title">Promotheus Nextflow metro map</title><desc id="desc">Validated inputs feed independent orthology and RNA-Seq routes. Merged orthology also feeds Microarray analysis. All routes converge at integration, followed by output validation. Rockhopper is optional within RNA_ANALYZE.</desc>', '<rect width="2000" height="820" fill="#282828"/>', '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#ededed}.title{font-size:38px;font-weight:700}.section{font-size:20px;font-weight:700;fill:#ccc}.label{font-size:20px;font-weight:700}.detail{font-size:16px;fill:#bdbdbd}.code{font-family:monospace;font-size:14px;fill:#aaa}</style>']
def text(x,y,t,cls='label',anchor='middle'):
 s.append(f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}">{escape(t)}</text>')
def panel(x,y,w,h,n,title):
 s.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="#353535" stroke="#515151"/>')
 text(x,y-15,f'{n}  {title}','section','start')
def path(d,c,width=5,dash=False):
 s.append(f'<path d="{d}" fill="none" stroke="{c}" stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"'+(' stroke-dasharray="8 7"' if dash else '')+'/>')
def station(x,y,title,code='',detail='',above=True):
 s.append(f'<rect x="{x-6}" y="{y-13}" width="12" height="26" rx="6" fill="#fff" stroke="#222" stroke-width="2"/>')
 text(x,y-31 if above else y+40,title)
 if code:text(x,y+43 if above else y+65,code,'code')
 if detail:text(x,y+69 if above else y+90,detail,'detail')
text(35,55,'Promotheus / Nextflow','title','start')
text(1960,53,'BACTERIAL TRANSCRIPTOMICS','detail','end')
panel(30,320,370,210,'1','Inputs')
panel(445,130,920,170,'2','Orthology · parallel by genome')
panel(445,350,920,150,'3','Microarray · parallel by GSE')
panel(445,550,920,170,'4','RNA-Seq · parallel download / alignment by GSE')
panel(1420,320,550,210,'5','Integration & validation')
green='#52c878';orange='#ff9a24';blue='#26aada';grey='#bdbdbd';purple='#c78cec'
# Validated input branches: RNA can proceed independently of orthology.
path('M85 415 H340',grey)
path('M340 415 V235 Q340 210 365 210 H1300',blue)
path('M340 415 V610 Q340 635 365 635 H1300',orange)
# Microarray receives validated shards and merged orthology.
path('M340 415 H755',grey,3,True)
path('M1170 210 H1310 Q1340 210 1340 240 V330 Q1340 340 1315 340 H730 Q705 340 705 365 V410 H1300',green)
# Orthology direct input to integration.
path('M1190 210 H1380 Q1400 210 1400 230 V389 Q1400 405 1420 405 H1490',blue)
path('M1300 410 H1490',green)
path('M1300 635 H1370 Q1390 635 1390 615 V435 Q1390 420 1410 420 H1490',orange)
# Merge all scientific routes before final verification.
path('M1490 405 H1530',blue)
path('M1490 410 H1530',green)
path('M1490 420 H1530',orange)
path('M1530 410 H1900',purple)
station(100,415,'XLSX × 2','','GEO datasets')
station(260,415,'Validate & split','VALIDATE_INPUTS')
station(575,210,'Genomes & proteomes','ORTHOLOGY_PREPARE','RefSeq / archived genomes')
station(870,210,'OrthoFinder','ORTHOLOGY_PAIR','Reference × each genome')
station(1170,210,'Merge gene IDs','ORTHOLOGY_MERGE','Wait for all comparisons')
station(930,410,'Chip & expression analysis','MICROARRAY_GSE','Requires merged orthology')
station(590,635,'Reads & reference','RNA_DOWNLOAD_GSE')
station(870,635,'Rsubread alignment','RNA_ALIGN_GSE')
station(1170,635,'Counts & expression','RNA_ANALYZE','featureCounts · DESeq2')
# Optional component within RNA_ANALYZE, not an independent scheduled process.
text(1170,586,'+ Rockhopper (optional)','detail')
station(1530,410,'Integrate','INTEGRATE_RESULTS')
station(1750,410,'Verify','VERIFY_OUTPUTS')
s.append('<path d="M1880 388 H1905 L1918 401 V438 H1880 Z" fill="#eee"/>')
text(1899,421,'TSV','code')
text(1899,471,'XLSX / TSV','detail')
# Legend and scope notes.
for x,c,t in [(45,blue,'Orthology'),(250,green,'Microarray'),(465,orange,'RNA-Seq'),(670,purple,'Final results')]:
 path(f'M{x} 770 h32',c);text(x+45,776,t,'detail','start')
text(1020,765,'Operons enabled: RNA_ANALYZE groups GSEs by organism.','detail','start')
text(1020,791,'Stations = Nextflow processes; dashed line = validated input shards.','detail','start')
s.append('</svg>')
Path('docs/images/promotheus-nextflow-dag.svg').write_text('\n'.join(s)+'\n')
