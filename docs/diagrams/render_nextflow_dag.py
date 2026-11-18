#!/usr/bin/env python3
"""Draw the scientific workflow with metro routes and reader-facing labels.

The layout follows the ten processes in main.nf. Reader-facing steps are accompanied by brief data-source and tool references;
process identifiers and configuration details remain in the technical guide.
"""
from html import escape
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
W,H=2000,1120
C=dict(orth='#008fcb',ma='#00ad62',rna='#f39a00',out='#9254d5',input='#66717e')
s=[f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="title description">
<title id="title">Promotheus workflow</title>
<desc id="description">Check selected studies and compare reference genomes to match genes. Download and analyze Microarray and RNA-seq data along separate routes. Gene matches support Microarray analysis and the final integration. Combine studies, check results and produce tables of gene expression, variability and promoter sequences. Optional operon analysis identifies genes transcribed together.</desc>
<rect width="2000" height="1120" fill="white"/>
<style>
text {{font-family:Arial,Helvetica,sans-serif;fill:#252525}}
.title {{font-size:30px;font-weight:700}}
.subtitle {{font-size:19px;fill:#555}}
.heading {{font-size:21px;font-weight:700}}
.name {{font-size:21px;font-weight:700}}
.body {{font-size:18px;fill:#454545}}
.note {{font-size:17px;fill:#666}}
.tech {{font-size:16px;fill:#666}}
</style>''']
def text(x,y,t,cls='body',anchor='middle',color=None):
    style=f' style="fill:{color}"' if color else ''
    s.append(f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}"{style}>{escape(t)}</text>')
def rect(x,y,w,h,stroke='#d2d2d2',fill='none'):
    s.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" stroke="{stroke}" fill="{fill}"/>')
def panel(x,y,w,h,n,label):
    rect(x,y,w,h);text(x+20,y+37,f'{n}. {label}','heading','start')
def route(d,color,dashed=False):
    dash=' stroke-dasharray="8 10"' if dashed else ''
    s.append(f'<path d="{d}" stroke="{color}" stroke-width="6" fill="none" stroke-linecap="round" stroke-linejoin="round"{dash}/>')
def arrow(x,y,color,down=False):
    points=f'{x-8},{y-10} {x},{y} {x+8},{y-10}' if down else f'{x-10},{y-8} {x},{y} {x-10},{y+8}'
    s.append(f'<polyline points="{points}" stroke="{color}" stroke-width="3.5" fill="none"/>')
def station(x,y,h=26):
    s.append(f'<rect x="{x-6}" y="{y-h/2}" width="12" height="{h}" rx="5" stroke="#58616b" stroke-width="1.5" fill="#fff"/>')
def node(x,y,title,detail):
    station(x,y);text(x,y-32,title,'name');text(x,y+44,detail)

text(40,65,'Promotheus workflow','title','start')
text(40,101,'From public gene expression studies to promoter profiles','subtitle','start')
panel(40,420,310,220,1,'Select studies')
panel(440,170,990,220,2,'Compare genomes')
text(1400,207,'Genomes can be compared in parallel','note','end')
panel(1030,460,400,180,3,'Microarray data')
panel(440,700,990,240,4,'RNA-seq data')
text(1400,737,'Studies can be processed in parallel','note','end')
panel(1530,170,430,770,5,'Combine & check results')

route('M90 540 H380',C['input'])
route('M380 540 V330 Q380 310 400 310 H1450 Q1470 310 1470 330',C['orth'])
route('M380 540 H1260',C['ma'],True)
route('M380 540 V800 Q380 820 400 820 H1450 Q1470 820 1470 800 V572 Q1470 552 1490 552 H1630',C['rna'])
route('M1260 310 V540',C['orth'])
route('M1470 330 V508 Q1470 528 1490 528 H1630',C['orth'])
route('M1260 540 H1630',C['ma'])
route('M1630 540 H1900 Q1920 540 1920 560 V672 Q1920 692 1900 692 H1760 Q1740 692 1740 712 V735',C['out'])
for x,y,k in [(415,310,'orth'),(760,310,'orth'),(1100,310,'orth'),(1260,433,'orth'),(1565,528,'orth'),(985,540,'ma'),(1500,540,'ma'),(415,820,'rna'),(760,820,'rna'),(1100,820,'rna'),(1565,552,'rna'),(1760,540,'out')]:
    arrow(x,y,C[k],(x,y)==(1260,433))
arrow(1740,730,C['out'],True)

rect(73,519,34,42,'#777','#fff')
for y in [532,540,548]:s.append(f'<path d="M80 {y} h20" stroke="#999"/>')
station(240,540)
text(240,508,'Check study details','name')
text(195,594,'Select studies and species','note')
text(195,619,'NCBI GEO · series IDs (GSE)','tech')
node(570,310,'Prepare genomes','Retrieve reference genomes')
text(570,378,'NCBI RefSeq','tech')
node(900,310,'Find matching genes','Identify corresponding genes')
text(900,378,'OrthoFinder · orthology','tech')
station(1260,310);text(1260,278,'Combine gene matches','name')
text(1290,365,'Gene matches','note','start',C['orth'])
station(1260,540)
text(1230,580,'Analyze expression','name')
text(1230,607,'Download and normalize array data','note')
text(1230,630,'GEO supplementary files · expression analysis','tech')
node(570,820,'Download study data','Retrieve sequencing data')
text(570,889,'GEO / GSE · SRA reads (SRR)','tech')
node(900,820,'Map reads to the genome','Assign reads to genomic locations')
text(900,889,'Rsubread','tech')
node(1260,820,'Measure gene expression','Estimate how active each gene is')
text(1260,889,'featureCounts · DESeq2','tech')
text(1260,919,'Rockhopper · optional operon detection¹','tech')

text(1560,263,'Bring studies together','name','start')
for y,t in [(298,'Link results through matching genes'),(329,'Correct study effects (ComBat)'),(360,'Summarize expression and variability'),(391,'Retrieve promoter sequences')]:
    s.append(f'<circle cx="1564" cy="{y-6}" r="2.5" fill="#666"/>');text(1577,y,t,'note','start')
station(1630,540,48);text(1630,493,'Combine','name')
station(1850,540);text(1850,493,'Check results','name')
text(1745,623,'Check completeness and consistency','note')
rect(1560,735,370,175,'#999','#fff')
text(1584,774,'Results','name','start')
text(1584,809,'Expression, percentiles and CV','body','start')
text(1584,840,'Promoter sequences','body','start')
text(1584,880,'XLSX / TSV + validation report','note','start')
s.append('<path d="M40 989 H1960" stroke="#bcbcbc" stroke-width="1"/>')
for x,k,label in [(40,'orth','Genome comparison'),(310,'ma','Microarray'),(520,'rna','RNA-seq'),(720,'out','Combined results')]:
    route(f'M{x} 1028 h30',C[k]);text(x+43,1034,label,'note','start')
station(1025,1028);text(1043,1034,'Analysis step','note','start')
route('M1260 1028 h34',C['ma'],True);text(1306,1034,'Selected studies','note','start')
text(40,1065,'GSE: GEO study accession.   SRA: sequencing-read archive.   SRR: SRA run accession.','tech','start')
text(40,1094,'¹ Operons are groups of genes transcribed together. Their detection is optional within RNA-seq analysis.','note','start')
s.append('</svg>')
(ROOT/'docs/images/promotheus-nextflow-dag.svg').write_text('\n'.join(s)+'\n')
