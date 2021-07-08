A list that contains information about the software maintainer(s). Might be convenient to add an
anchor (with '&' symbol) to this field if the maintainer(s) is also the also the tool creator so they can also be 
referenced in the optional `creator` field. `identifier` fields are not required, but if included should be a unique identifier such as orcid.
~~~yaml
contactPoint: 
  - &Jane
    name: Jane Schmoe 
    email: jane@theschmoes.com
    identifier: https://orcid.org/0000-0001-6022-9825
  - name: Joe Schmoe
    email: joe@theschmoes.com
~~~