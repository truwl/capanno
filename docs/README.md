# Docs

Markdown files that include component files and final docs that are assembled with templates (in `/templates` directory. To produce assembled docs
use [pandoc](https://pandoc.org) with the  [pandoc-include](https://pypi.org/project/pandoc-include/) filter. e.g.

~~~bash
$ pandoc docs/templates/CommandLineTool_metadata_guide.md --filter pandoc-include -o docs/CommandLineTool_metadata.md -t markdown -s
~~~