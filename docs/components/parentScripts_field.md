List of scripts that the script uses (imports of other scripts). Used to create navigable relationships.

ex.

```yaml
parentScripts:
  - name: 
    alternateName:  # useful for looking up scripts if searching 'name' does not provide hits.
    softwareVersion:  # Needed to find correct version of script.
    identifier: # truwl identifer if it exists. If this is provided, other fields not used.
```