Relative path (or list of paths) to metadata file(s) that contains fields that the primary metadata should inherit from. Fields from parentMetadata are only used if values are not provided in the primary file and metaadata listed earlier takes precedence over metadata listed later.

ex.

```yaml
parentMetadata: ../common/foo-metadata.yaml
```

or 

```yaml
parentMetadata:
  - ../common/foo-metadata.yaml  # fields in foo-metadata.yaml take precedence to fields in bar-metadata.yaml
  - ../common/bar-metadata.yaml
```