# Docs

Directory for markdown files that include component files, template files, and final assembles doc files. To produce assembled docs from template
use [pandoc](https://pandoc.org) with the  [pandoc-include](https://pypi.org/project/pandoc-include/) filter. e.g.

~~~bash
$ pandoc docs/templates/Tools.md --filter pandoc-include -o docs/ToolGuide.md -t markdown -s
~~~