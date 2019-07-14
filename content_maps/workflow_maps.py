from collections import OrderedDict
from itertools import chain

example_workflows = OrderedDict([
    ('WF_1f4d8f.cb', 'cwl-workflows/example_workflows/cat_sort/1.0/cat_sort.cwl'),
    ('WF_0b650e.e4', 'cwl-workflows/GA4GH_Workflow_Execution_Challenge/hello_world/1.0/hello_world_with_checker_workflow.cwl'),
    ('WF_680a59.6e', 'cwl-workflows/GA4GH_Workflow_Execution_Challenge/md5sum_workflow/master__42b7f70/md5sum-workflow.cwl'),
    ('WF_ba9852.6e', 'cwl-workflows/GA4GH_Workflow_Execution_Challenge/md5sum_checker/master__42b7f70/checker-workflow-wrapping-workflow.cwl')
])


everything = OrderedDict(chain(example_workflows.items()))