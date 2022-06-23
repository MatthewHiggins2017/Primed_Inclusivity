# **Primed Inclusivity**
-----------------------

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/69949ec20c6a40d5be32ae5fe00115ee)](https://www.codacy.com/gh/MatthewHiggins2017/Primed_Inclusivity/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=MatthewHiggins2017/Primed_Inclusivity&amp;utm_campaign=Badge_Grade)
-----------------------

**Note**-*This page is still under development*

To Do:

* Provide detailed breakdown of Output Guide
* Provide detailed breakdown of Representation Guide

-----------------------
**Summary**


Primed Inclusivity is a python-based programmable framework, enabling researchers to detect and gauge the impact of primer binding site heterogeneity on NAAT-based assay performance.

-----------------------
**Installation**


From GitHub.

```
git clone https://github.com/MatthewHiggins2017/Primed_Inclusivity.git
cd ./Primed_Inclusivity/
conda env create -f PIEnv.yml
conda activate PIEnv
python setup.py install

```


From Conda

```
# Coming Soon
```

-----------------------
**Example Analysis**

Please navigate to the wiki where you will find a tutorial walk through of how to use the PrimedInclusivity tool and add custom modules. [**Click Here.**](https://github.com/MatthewHiggins2017/Primed_Inclusivity/wiki/Tutorial)  

------------------------

**Input Parameters**


| **Variable ID**                   | **Symbol**  | **Description** |
|--------------------------------   |-----    |-------------|
| RunID                             | -i      |  Identification ID for files generated.                                                     |
| Threads                           | -@      |  Number of simultaneous processes to run.                                                   |
| VCF                               | -v      |  Path to input multisample  VCF representing the population                                 |
| SamplesFile                       | -SF     |  NAN                                                                                        |
| Target                            | -t      |  Path to input fasta file of target organisms which primers are supposed to bind.           |
| Primers                           | -p      |  Path to primer input file. See section below for more information.                         |
| Sets                              | -s      |  Path to sets input file. See section below for more information.                           |
| Output                            | -o      |  Path to output directory. directory                                                                   |
| RepresentationGuide               | -RG     |  Path to representation guide. NAN                                                                                        |
| OutputGuide                       | -OG     | Path to output guide NAN                                                                                        |           
| OutputType                        | -OT     | Define what information should be written from database to CSV files. Options: All (A), Sets (S), Binding Classification (BC), Sample vs Primer (SvP), Representation Tables (RT). Ensure input values are comma seperated.                                                                                      |           
| DatabasePath                      | -DB     | Path to Previously Generated Database. NAN                                                                                        |          
| UpdateProbability                 | -UP     |  Using existing binding sites update output values accordingly.                                                                                         |           
| AnchorRegion                      | NA      |  See Optimising Primer Binding. Github Wiki.                                                                                         |           
| AnchorRegionMismatchScaler        | NA      |  See Optimising Primer Binding. Github Wiki.                                                                                         |           
| AnchorAdjacentMismatchScalar      | NA      |  See Optimising Primer Binding. Github Wiki.                                                                                         |           
| AdjacentMismatchScalar            | NA      |  See Optimising Primer Binding. Github Wiki.                                                                                         |          
| ComplementaryScalar               | NA      |  See Optimising Primer Binding. Github Wiki.                                                                                         |          
| MismatchScaler                    | NA      |  See Optimising Primer Binding. Github Wiki.                                                                                         |                    


-------------------------------
