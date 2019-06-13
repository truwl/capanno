from collections import OrderedDict
from itertools import chain

example_workflows = OrderedDict([
    ('WF_1f4d8f.cb', 'workflows/example_workflows/cat_sort/1.0/cat_sort.cwl'),
    ('WF_0b650e.e4', 'workflows/example_workflows/hello_world/1/hello_world_with_checker_workflow.cwl'),
    ('WF_680a59.6e', 'workflows/example_workflows/md5sum_workflow/1/md5sum-workflow.cwl'),
    ('WF_ba9852.6e', 'workflows/example_workflows/md5sum_checker/1/checker-workflow-wrapping-workflow.cwl')
])


everything = OrderedDict(chain(example_workflows.items()))