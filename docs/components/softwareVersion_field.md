The `softwareVersion` field has two subfields: `versionName` and `indluededVersions`. `versionName` should follow the tool's version convention as much as possible and can include variable portions, e.g. v3.x.y. For some version conventions it may be appropriate for the versionName to be a range, e.g. v301-v399 

Example:
```yaml
softwareVersion: 
  versionName: 2.x
  includedVersions:
    - 2.1.1
    - 2.1.2
    - 2.2.1
    - 2.2.2
```