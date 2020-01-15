# Understand cancer evolution through single cell expression dynamics

Sehyun Oh (sehyun.oh@sph.cuny.edu) 
<br/>
Ludwig Geistlinger (ludwig.geistlinger@sph.cuny.edu) 
<br/>
Ola Oni (olo4002@med.cornell.edu)

Department of Epidemiology and Biostatistics    
CUNY School of Public Health   
New York, NY 10027

Nur-Taz Rahman (nur-taz.rahman@yale.edu)
<br/>
Simbonis Bioinformatics Fellow
Cushing/Whitney Medical Library
333 Cedar St, New Haven, CT 06511

Alejandro Jimenez-Sanchez (ajs.scientia@gmail.com)
<br/>
Postdoctoral Research Fellow
The Dana Pe'er Lab
Computational & Systems Biology Program
Sloan Kettering Institute, NY
New York

## Question
Can we detect the direction of cancer evolution by analyzing the expression dynamics? If cancer samples are a mixture of different cells (cancer subtypes and/or tissues of origin), can we infer the stage and/or origin of cells from their evolutionary trajectory?


## Approach
#### Tasks
- Cell type identification    
- Known cancer subtypes or tissues of origin   
- Identification of subclonal cell populations: using mutations in RNA-seq data rather than just gene expression?     
- Evolutionary order of subclones and observation of EMT by trajectory analysis   
- Study a priori markers of metastatic potential and whether they are early or late events   

#### Deliverable
- Analysis workflows for each steps → make it reusable for new datasets
- Compare trajectory inference tools on different cancer samples → Benchmark summary


## Tools
#### Background
- [Psudotime Cell Trajectories](https://docs.google.com/presentation/d/e/2PACX-1vQuzaq2kbvEEc3mrUwILcCHuovrKKZWU45EQVEzWISgRVgl3A5KYR1FuY1cS2w0DHG-0wO19zGtvaNj/embed?start=false&loop=false&delayms=3000&slide=id.p) (a part of Broad Institute workshop material)

- List of Trajectory Inference Methods (https://github.com/dynverse/dynmethods#list-of-included-methods) and Benchmarking package (https://dynverse.org/)

- Functional Pseudotime analysis (https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html), originally from [the Hemberg lab course material](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#pseudotime-analysis)

- [A comparison of single-cell trajectory inference methods](https://www.nature.com/articles/s41587-019-0071-9)

#### For codeathon
1. [slingshot](https://bioconductor.org/packages/release/bioc/html/slingshot.html)
2. [monocle3](https://cole-trapnell-lab.github.io/monocle3/)

## Datasets
1. _**[tissue of origin]**_ Identification of grade and origin specific cell populations in serous epithelial ovarian cancer by single cell RNA-seq
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118828) [[GSE118828](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118828)]

2. _**[primary vs. metastasis]**_ Single-Cell Transcriptomic Analysis of Primary and Metastatic Tumor Ecosystems in Head and Neck Cancer
(https://www.cell.com/cell/fulltext/S0092-8674(17)31270-9?cid=tw%26p) [[GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322)]

3. _**[longitudinal]**_ Longitudinal single-cell RNA sequencing of patient-derived primary cells reveals drug-induced infidelity in stem cell hierarchy
(https://www.nature.com/articles/s41467-018-07261-3) [[GSE117872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117872)]


## Misc.
#### Other tools  
Below are the examples of trajectory inference tools on cancer, with very poor documentation - not enough to reproduce. 

- Mapping lung cancer epithelial-mesenchymal transition states and trajectories with single-cell resolution
(https://www.nature.com/articles/s41467-019-13441-6) --> TRACER

- High Resolution Comparison of Cancer-Related Developmental Processes Using Trajectory Alignment 
(https://www.biorxiv.org/content/10.1101/469601v1.full) --> devMap

