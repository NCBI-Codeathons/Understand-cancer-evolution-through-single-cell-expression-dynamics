# Understand cancer evolution through single cell expression dynamics
## Team
**Sehyun Oh** (sehyun.oh@sph.cuny.edu, Lead)
<br/>
**Ludwig Geistlinger **(ludwig.geistlinger@sph.cuny.edu)

Department of Epidemiology and Biostatistics    
CUNY School of Public Health   
New York, NY 10027

**Nur-Taz Rahman** (nur-taz.rahman@yale.edu)
<br/>
Simbonis Bioinformatics Fellow
<br>
Cushing/Whitney Medical Library
<br>
333 Cedar St, New Haven, CT 06511

**Alejandro Jimenez-Sanchez** (ajs.scientia@gmail.com)<br/>
Postdoctoral Research Fellow<br/>
The Dana Pe'er Lab<br/>
Computational & Systems Biology Program<br/>
Sloan Kettering Institute, NY<br/>
New York, NY 10065<br/>

**Ola Oni** (olo4002@med.cornell.edu)
<br/>
Tri-I PhD Program in Computational Biology and Medicine
<br/>
Weill Cornell Medicine
<br/>
New York, NY 10021

## Question
Can we understand cancer stage using SC trajectory inference tools?

## Deliverables
- **Benchmark summary**: apply trajectory inference tools on different single-cell cancer samples    
- **Analysis workflows**    
- **Proposal on tool updates** (e.g. what is the biological difference between differentiating cells and cancer cells need to be considered)

## Datasets
R and Pyhton scripts for downloading these datasets are available under `datasets` folder.   

1. _**[tissue of origin]**_ Identification of grade and origin specific cell populations in serous epithelial ovarian cancer by single cell RNA-seq
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118828) [[GSE118828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118828)]

2. _**[primary vs. metastasis]**_ Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer
(https://www.cell.com/cell/fulltext/S0092-8674(17)31270-9?cid=tw%26p) [[GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322)]


## Tools
- [slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html)
- [monocle3](https://cole-trapnell-lab.github.io/monocle3/)
- [Plantir](https://github.com/dpeerlab/Palantir)
- [PAGA](https://github.com/theislab/paga)
- [dyno](https://dynverse.org/)


## Lessons learned

## Next steps


## Misc.
#### Background
- [Psudotime Cell Trajectories](https://docs.google.com/presentation/d/e/2PACX-1vQuzaq2kbvEEc3mrUwILcCHuovrKKZWU45EQVEzWISgRVgl3A5KYR1FuY1cS2w0DHG-0wO19zGtvaNj/embed?start=false&loop=false&delayms=3000&slide=id.p) (a part of Broad Institute workshop material)

- List of Trajectory Inference Methods (https://github.com/dynverse/dynmethods#list-of-included-methods) and Benchmarking package (https://dynverse.org/)

- Functional Pseudotime analysis (https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html), originally from [the Hemberg lab course material](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#pseudotime-analysis)

- [A comparison of single-cell trajectory inference methods](https://www.nature.com/articles/s41587-019-0071-9)

#### Other tools  
Below are the examples of trajectory inference tools on cancer, with very poor documentation - not enough to reproduce.

- Mapping lung cancer epithelial-mesenchymal transition states and trajectories with single-cell resolution
(https://www.nature.com/articles/s41467-019-13441-6) --> TRACER

- High Resolution Comparison of Cancer-Related Developmental Processes Using Trajectory Alignment
(https://www.biorxiv.org/content/10.1101/469601v1.full) --> devMap
