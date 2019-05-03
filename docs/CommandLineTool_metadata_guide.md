Metadata Guide
==============
### Introduction

This documents how to describe metadata for command line tools described by the common workflow language (CWL) in this repository. 
Separate metadata files must be specified for each version of a tool (excluding patch versions). Metadata files are
specified in YAML. Most metadata keys are defined by [schema.org](https://schema.org/) vocabularies and can be adapted to be 
imported by CWL files.
### Metafile types 
There are three types of metadata files:
-  **[Complete tool metadata](#complete):** 
Metadata file for a tool that is not separated into separate subtools. 
-  **[Parent metadata](#parent)** 
 Metadata that is inherited by multiple subtools.
- **[Subtool metadata](#subtools)**
Metadata that is specific to a subtool. 

Templates are provided for each of these types of metadata in the [templates](../templates/) directory.

### List order
Metadata fields that accept lists as values should be ordered by importance if applicable. 
### Metadata sources
We recommend populating metadata files from [bio.tools](https://bio.tools/) if available. A python script to pre-populate 
metadata files from the bio.tools api should be coming soon and will be made available.

Metadata can also be obtained from tool manual and help pages, [SciCrunch](https://scicrunch.org/), or other web resources.


## <a name="complete"><a/>Complete tool metadata
Metadata file fields for a tool that is not separated into subtools.

### Required Fields
#### <a name="name1"><a/>name
The name of the tool. Typically, extensions are left out of the name. e.g. 'Picard', not 'picard.jar' although this is
not a requirement. Alternate names can optionally be stored in the [alternateName](#alternatename) list. 
ex: `name: Picard`
#### <a name="sversion1"></a>softwareVersion
A string that specifies the version of the tool that the CWL file is valid for. Versions must follow the 
version conventions used by the tool.
ex: `softwareVersion: v3.2`

#### <a name="version1"><a/>version
Specifies the version of the CWL file. 
Must follow [semantic versioning](https://semver.org/spec/v2.0.0.html) conventions.
A 1.0 or greater version of a CWL document must contain valid CWL syntax, 
follow the required [best practices](CommandLineTool_guide.md), and 
describe enough input parameters to be useful to others. Breaking changes that constitute a major revision
include changing the name/id of a command input parameter that would make any job files not run properly.
ex: `version: 0.2.1`


### Recommended Fields
#### <a name="desc1"><a/>description
Description of the tool. Should be taken directly from the tool documentation if it exists.

ex:
~~~yaml
description: |
   grep  searches  for  PATTERN  in  each  FILE.  A FILE of “-” stands for standard input.  If no FILE is given, 
   recursive searches examine the working directory, and nonrecursive searches read standard input.  By default, 
   grep prints the matching lines.

   In addition, the variant programs egrep, fgrep and rgrep are the same as grep -E, grep -F, and grep -r, respectively.  
   These variants are deprecated, but are provided for backward compatibility.
~~~
#### <a name="repo1"><a/>codeRepository
Code repository information. `URL` should specify a web page (should not have a .git extension).

ex: cwltool
~~~yaml
codeRepository:
  name: GitHub
  URL: https://github.com/common-workflow-language/cwltool

~~~
#### <a name="license1"><a/>license
Software license for the tool defined by [SPDX identifier](https://spdx.org/licenses/). e.g. `license: Apache-2.0`
#### <a name="site1"><a/>Website
A list of websites associated with the tool, excluding the code repository which should be provided in `codeRepository`.

ex:
~~~yaml
Website:
  - name: Samtools  # The name of the website.
    description: Samtools homepage.
    URL: http://www.htslib.org/
  - name: Wikipeda
    description: Wikipedia entry for SAMtools
    URL: https://en.wikipedia.org/wiki/SAMtools

~~~

#### <a name="cpoint1"><a/>contactPoint
A list that contains information about the software maintainer(s). Might be convenient to add an
anchor (with '&' symbol) to this field if the maintainer(s) is also the also the tool creator so they can also be 
referenced in the optional `creator` field. `identifier` fields are not required, but recommended if available.

~~~yaml
contactPoint: 
  - &Jane
    name: Jane Schmoe 
    email: jane@theschmoes.com
    identifier: https://orcid.org/0000-0001-6022-9825
  - name: Joe Schmoe
    email: joe@theschmoes.com
~~~
 
#### <a name="pub1"><a/>publication
A list that describes publications related to the tool. The first entry should be the reference to cite the tool. 

ex:
~~~yaml
publication:
  - headline: 'BEDTools: a flexible suite of utilities for comparing genomic features.' # Title goes here.
    identifier: 10.1093/bioinformatics/btq033  # DOI goes here
  - headline: 'BEDTools: The Swiss-Army Tool for Genome Feature Analysis.'
    identifier: 10.1002/0471250953.bi1112s47
~~~
#### <a name="key1"></a>keywords
List of tags to categorize the tool by topic and operation specified with an 
[edam](http://bioportal.bioontology.org/ontologies/EDAM?p=classes)  or other ontology identifier.
EDAM keywords are preferred. If you wish to provide a keyword that is not in an ontology, it may be 
specified with the keys `name` and `category`. The value of `category` must be either 'topic' or 
'operation'. These categories have the meanings defined by [EDAM](http://edamontology.org/page).

ex.
~~~yaml
keywords:
  - http://edamontology.org/operation_3182  # Genome alignment
  - http://edamontology.org/topic_0085  # Functional Genomics
  - name: ENCODE
    category: topic
~~~

### Optional Fields
#### <a name="altname1"><a/>alternateName
List of alternate names for the tool. This is convenient place to put any other names for the tool that someone 
might use to search for it.
ex: `alternateName: []`
#### <a name="creator1"><a/>creator
List of the tool's creator(s). Might be redundant if this is captured in `contactPoint` and/or `publication`.
Can alias `contactPoint` if this is the case.
~~~yaml
creator:
  - *Jane # alias to &Jane in the contactPoint field.
  - name: Bob Bobbins
    email: bob@bobsbobbins.eu
    identifier: 
~~~
#### <a name="plang1"><a/>programmingLanguage
List of programming languages that the tool is written in `programmingLanguage: [Python2, C]`
#### <a name="datepub1"><a/>datePublished
The release date of the particular version of software in international standard date notation (YYYY-MM-DD). 
ex: `datePublished: 2003-04-14`
#### <a name="dload1"><a/>downloadURL
URL to download the tool. `downloadURL: https://github.com/arq5x/bedtools2/archive/v2.27.1.zip`


## <a name="parent"><a/>Parent (common) Metadata
Metadata file fields for metadata that is inherited by multiple subtools.

### Required Fields
#### name
Name of the main tool. ex: `pip`. Also see [name](#name1) above.
#### softwareVersion
Same as [softwareVersion](#sversion1) above.
#### featureList
Required for tools that are divided into subtools. The names in the list must correspond to `subtoolName` 
provided in the metadata file for the subtool, if one has been created, and should be the string used to identify the 
subtool when running it from the command line, excluding any dashes. `featureList` must include all the subtools of the tool
even when a metadata file or CWL file has not been created for the subtool. When the `featureList` field is populated 
the metadata file will not correspond to a single CWL file but will be used as a common metadata file for multiple 
subtools to inherit from and must be placed in the tool's `common/` directory.

ex. For pip:
```yaml
featureList:
  - install
  - download
  - uninstall
  - freeze
  - list
  - show
  - check
  - config
  - search
  - wheel
  - hash
  - completion
  - help  # This one is not really necessary
```

### Recommended Fields
#### description
Description of the tool that applies to the tool as a whole, not to an individual subtool. 
Also see [description](#desc1) above.
#### codeRepository
Same as [codeRepository](#repo1) above.
#### license
Same as [license](#license1) above.
#### Website
Same as [Website](#site1) above.
#### contactPoint
Same as [contactPoint](#cpoint1) above.
#### publication
Same as [publication](#pub1) above.
#### keywords
Terms that are relevant to all subtools of the tool. See [keywords](#key1) above.

### Optional Fields
#### alternateName
Same as [alternateName](#altname) above.
#### creator
Same as [creator](#creator1) above.
#### programmingLanguage
Same as [programmingLanguage](#plang1) above.
#### datePublished
Same as [datePublished](#datepub1) above.
#### downloadURL
Same as [downloadURL](#dload1) above.



## <a name="subtools"><a/>Subtool metadata
Metadata file fields for a CWL file that describes a subtool of a tool.
### Required Fields
#### applicationSuite
Information to identify the main tool. The subtool will inherit metadata from the parent metadata.
The `name` and `SoftwareVersion` fields must match the corresponding fields in the [main tool metadata file].(#)
~~~yaml
applicationSuite:
  name: pip
  softwareVersion: v19.0
  identifier:  # truwl identifier of main tool metadata, if it exists.
~~~
#### name
This must correspond to the name of the subtool as specified in the [featureList](#featurelist) field of the primary tool metadata file.
ex: `search`
#### version
Same as [version](#version1) above. Each subtool CWL file
must have its own version.
### Recommended Fields
#### description
Description of the subtool. Should contain subtool specific information.

ex. pip search
~~~yaml
description:  |
  Search PyPI for packages.
~~~
##### keywords
Terms that are specfic to the subtool. 
Terms that are relevant to all subtools of a tool should be defined in the common metadata file.
See [keywords](#keywords) above.

### Optional Fields
#### alternateName
List of alternate names for the subtool. 
ex: gzip -c . gzip is parent tool name. 'c' is subtool name. zcat is an alternateName.
~~~yaml
alternateName:
  - zcat
~~~
