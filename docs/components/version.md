#### <a name="version1"><a/>version

Specifies the version of the CWL file. 
Must follow [semantic versioning](https://semver.org/spec/v2.0.0.html) conventions.
A 1.0 or greater version of a CWL document must contain valid CWL syntax, 
follow the required [best practices](CommandLineTool_guide.md), and 
describe enough input parameters to be useful to others. Breaking changes that constitute a major revision
include changing the name/id of a command input parameter that would make any job files not run properly.
ex: `version: 0.2.1`