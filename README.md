# **PrimedInclusivity**
-----------------------

**Note**-*This page is still under development*

To Do:

* Create quick start example
* Provide detailed breakdown of Classification Engine
* Provide detailed breakdown of Output Engine
* Provide detailed breakdown of Output Guide
* Provide detailed breakdown of Representation Guide

-----------------------
### **Summary**


PrimedInclusivity is a python-based programmable framework, enabling researchers to detect and gauge the impact of primer binding site heterogeneity on NAAT-based assay performance.




-----------------------
### **Installation**


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


------------------------

# **Input Parameters**


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
| UpdateProbability                 | -UP     |  NAN                                                                                        |           
| AnchorRegion                      | NA      |  NAN                                                                                        |           
| AnchorRegionMismatchScaler        | NA      |  NAN                                                                                        |           
| AnchorAdjacentMismatchScalar      | NA      |  NAN                                                                                        |           
| AdjacentMismatchScalar            | NA      |  NAN                                                                                        |          
| ComplementaryScalar               | NA      |  NAN                                                                                        |          
| MismatchScaler                    | NA      |  NAN                                                                                        |           
| ThermoMinLength                   | NA      |  NAN                                                                                        |           





-------------------------------
