## EmbryoCV is an open-source Python package for the analysis of video describing embryonic development

EmbryoCV is a Python package developed by a team of scientists to quantify different biological aspects of developing organisms from video datasets. The package EmbryoCV was developed to analyse the large long-term image datasets acquired using an open-source video microscope OpenVIM www.openvim.org. OpenVIM can generate high quality long-term video datasets of large numbers of aquatic embryos for extended periods of time. Consequently OpenVIM is a powerful tool for visualising and quanitfying the dynamic period of embryonic development in a way that is simply not possible with manual observation. 

However, the downside of this capability is the generation of vast quantities of video data and the extremely time consuming nature of manual analysis. It therefore becomes necessary to quanitfy only particular aspects of biological development, or to limit the scale of experiments for methodological, rather than biological reasons. EmbryoCV is the solution to this issue. It provides a high-throughput pipeline for measuring phenome-level responses of aquatic embryos with unprecedented temporal, spatial and functional resolution.

### What can EmbryoCV do?
Biological development is complex and dynamic. EmbryoCV has been designed to capture as much of this complexity as possible in the form of biologically relevant measures. These measurements form a powerful phenome-level dataset and owing to the capacity of OpenVIM to record many hundred embryos simultaneously, whilst simultaneously controlling the embryonic environment, these responses can be quantified in large numbers of aquatic embryos and in different environmental conditions.

[Scrolling data gif] (assets/embryocvScrollingData.gif)

<img src="assets/embryocvScrollingData.gif" width="200" height="100" />
### How do I use EmbryoCV?
EmbryoCV is intentionally simple to use, consisting of a small number of user callable functions. A user begins by setting up an 'instance' of EmbryoCV for their analysis. They must provide EmbryoCV the ...
..
..


### How is EmbryoCV structured?
EmbryoCV has a modular internal structure, consisting of the modules dataHandling, imageAnalysis, eggUI, dataIntegration, dataAnalysis and the parent module EmbryoCV. However, this complexity is irrelevant to the user [More detailed information on the structure of EmbryoCV](programStructure.md)



