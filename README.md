# CCA
* Canonical Correlates Analysis (CCA) is a multivariate analysis technique for maximizing the correlations between orthagonalized sets of independent and dependent variates
* The codes provided here are for the CCA as performed in Perry et al., (2017, *in review*), where
  + CCA linked rs-fMRI patterns to demographic and cognitive measures in [101 cognitively-normal elders](https://cheba.unsw.edu.au/project/sydney-memory-and-ageing-study)
* And are modified from [Smith et al's., (2015)](http://www.nature.com/neuro/journal/v18/n11/full/nn.4125.html) HCP investigation 

![edgebund](https://cloud.githubusercontent.com/assets/23748735/26387350/3a32b8f8-4090-11e7-967a-ac9ba5a72c7a.png)

*Adapted from Perry et al., (2017)

### Codes included within repository: 
* Normalization and decomposition of functional networks
* CCA
* Basic visualization output
* Parcellation templates used in functional network construction

### Required dependencies
* MatLab
* [PALM](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM)
* [FSLNETS](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets)
* [BrainNetViewer](https://www.nitrc.org/projects/bnv/) (optional for visualisation)

## Getting started
### 1. Data required:
* Functional network matrices of all subjects (i.e. dependent variates in CCA) (*connectivity matrices*)
* Design Matrix of non-imaging measures (i.e. independent variates in CCA) (*DM*)
* Motion parameters (i.e. framewise displacement) of functional images (*motionFD*)
* Centroids of parcellation template employed in functional network construction (*COG*)

### 2. Performing the CCA
* Within a MatLab terminal run:
  + `[CCAout] = cca_functional(connectivitymatrices, DM, motionFD, COG)`

### 3. Extracting CCA results
* The resultant data will be stored within the Matlab structure `CCAout`
* Which stores important information in the fields of `CCAout`, such as:
  + *grotU*: Individual subjects weights for non-imaging measures captured by each CCA mode 
  + *grotV*: Individual subjects weights for functional connectivity patterns captured by each CCA mode
  + *grotR*: Correlation between the orthagonalized non-imaging and connectivity patterns
  + *conload*: Loadings of non-imaging measures onto each modes connectivity patterns (i.e. grotV)
  + *grotstats*: Parametric statistical output

### 4. Visualising the CCA output
* Extracted are the connectivity edges and nodes that are most strongly expressed (i.e. top 250 connections) by the first CCA mode
  + For both positive and negative associations with CCA mode
  + Users may want to modify their code, depending on their number of significant CCA modes
* Data is extracted in the `.nodes` and `.edge` format required for BrainNet Viewer:
* For example, for the top positive associations:
  + Nodes : `CCA_nodes250topposcons_mode1.nodes`
  + Edges : `CCA_250topposcons_mode1.edge`
  
![bnvhelpgithub](https://cloud.githubusercontent.com/assets/23748735/26387399/b53879e8-4090-11e7-977a-de4ea478fab8.png)

![hbm-mode1-poscons-axview](https://cloud.githubusercontent.com/assets/23748735/26387159/aafb3062-408e-11e7-9807-b862f3413ac1.png)  
For any questions and more advanced codes/data please contact Alistair Perry (QIMR Berghofer) (alistairgperry *at* gmail.com)
