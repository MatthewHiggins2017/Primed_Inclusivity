# **PrimedInclusivity**
-----------------------

**Note**-*This page is still under development*


-----------------------
## **Summary**
------------------------
PrimedInclusivity allows users to identify variants within a population which may affect primer binding and subsequently how they affect the probability of reaction success, exposing risk to false negatives.

PrimedInclusivity was developed due to the growing application of nucleic acid amplification tests (NAATs) in diagnostic and research settings. This includes Polymerase Chain Reaction (PCR), Loop Mediated Isothermal Amplification (LAMP), Recombinase Polymerase Amplification (RPA), etc.

As variant data is rarely homogenous across a given population, PrimedInclusivity has  been design to allow for the correction of any representation bias introduced from the availability of next-generation sequence data.

------------------------

# **Input Parameters**

| **Variable ID**                   | **Symbol**  | **Description** |
|--------------------------------   |-----    |-------------|
| RunID                             | -i      |  Identification ID for files generated.          |
| Threads                           | -@      |  Number of simultaneous processes to run.          |
| VCF                               | -v      |  Path to input multisample  VCF representing the population            |
| Target                            | -t      |  Path to input fasta file of target organisms which primers are supposed to bind.           |
| Primers                           | -p      |  Path to primer input file. See section below for more information.          |
| Sets                              | -s      |  Path to sets input file. See section below for more information.            |
| Output                            | -o      |  Path to output location.            |
| RepresentationJson                | -RS     |            |
| Haploid                           | -H      |             |
| AnchorLength                      | -AL     |             |
| AnchorRegionMismatchScaler        | -ARMS   |             |
| AnchorAdjacentMismatchScalar      | -ARAMS  |             |
| AdjacentMismatchScalar            | -AMS    |             |
| ComplementaryScalar               | -CS     |             |
| MismatchScaler                    | -MS     |             |
| Thermo                            | -T      |             |



-------------------------------


# **Representation Review**


The Representation json input file is split into two key dictionaries outlined below:

* RepTags - These are the user defined classification variables and corresponding marginal distribution for each relative category.

* SampleTags - These are classification assignments for each sample present. According to each category.
