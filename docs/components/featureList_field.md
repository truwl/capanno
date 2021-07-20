A list that contains the names of the subtools for a tool. The names in the list must correspond to `subtoolName` 
provided in the metadata file for the subtool, if one has been created, and should be the string used to identify the 
subtool when running it from the command line, excluding any dashes. If the tool can be called without a subcommand, the name of the main tool is `__main__` and should be included in the list.

e.g. for pip:
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