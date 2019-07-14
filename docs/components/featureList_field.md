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