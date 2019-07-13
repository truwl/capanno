List of tags to categorize the method by topic and operation specified with an 
[edam](http://bioportal.bioontology.org/ontologies/EDAM?p=classes)  or other ontology identifier.
EDAM keywords are preferred. If you wish to provide a keyword that is not in an ontology, it may be 
specified with the keys `name` and `category`. The value of `category` must be either 'topic' or 
'operation'. These categories have the meanings defined by [EDAM](http://edamontology.org/page).

ex.

~~~yaml
keywords:
  - http://edamontology.org/operation_3182  # Genome alignment
  - http://edamontology.org/topic_0085  # Functional Genomics
  - name: ENCODE
    category: topic
~~~
