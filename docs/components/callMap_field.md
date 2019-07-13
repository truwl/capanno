The `callMap` field relates the steps in a workflow to the tool, script, or or workflow that is called. All steps must be in the callMap list.

ex:

````yaml
callMap:
  - id:  # id (name) of the workflow step
    identifier: # truwl identifier of the tool, script, or workflow that is called.

````