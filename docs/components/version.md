Specifies the version of the CWL file and is used to determine its development state.
Must follow [semantic versioning](https://semver.org/spec/v2.0.0.html) conventions.
A 1.0 or greater version of a CWL document must contain valid CWL syntax, 
follow the required [best practices](CommandLineTool_guide.md), and 
describe enough input parameters to be useful to others. Breaking changes that constitute a major revision
include changing the name/id of a command input parameter that would make any job files not run properly. If not specified, will be initialized to 0.1.0
ex: `version: 0.2.1`