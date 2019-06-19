#### applicationSuite
Information to identify the main tool. The subtool will inherit metadata from the parent metadata.
The `name` and `SoftwareVersion` fields must match the corresponding fields in the [main tool metadata file].(#)

~~~yaml
applicationSuite:
  name: pip
  softwareVersion: v19.0
  identifier:  # truwl identifier of main tool metadata, if it exists.
~~~