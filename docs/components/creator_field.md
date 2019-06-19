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