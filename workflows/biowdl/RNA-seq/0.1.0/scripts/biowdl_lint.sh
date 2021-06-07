#!/usr/bin/env bash

# Validate every WDL file in the repository with womtool validate and miniwdl check.
# Also check parameter_meta presence unless the argument "skip-wdl-aid" is given.
set -e
for WDL_FILE in $(git ls-files *.wdl)
  do
    echo $WDL_FILE
    womtool validate $WDL_FILE
    miniwdl check $WDL_FILE
    grep Copyright $WDL_FILE || bash -c "echo No copyright header in $WDL_FILE && exit 1"
    # Run WDL-AID in strict mode (error if parameter_meta is missing for any
    # inputs). WDL-AID also errors if there is no workflow in the WDL file
    # but in this case we don't care about that. As such if WDL-AID errors
    # we check if it is the error we care about.
    wdl-aid --strict $WDL_FILE > /dev/null 2> wdl-aid_stderr ||
    if grep -z "ValueError: Missing parameter_meta for inputs:" wdl-aid_stderr
      then
        exit 1
    fi
  done

# For each submodule
# 1. Check if there is a tag. If so, we are using a stable version.
# 2. If not, we must be on the main develop branch. Make sure it is the
#    latest version of this branch.
# Afterwards check if there are changes. If so, some submodules have been
# updated to a newer commit, so that means they were not on develop. In that
# case exit 1.

git submodule foreach \
bash -c '\
if [ "$(git tag --contains)" == "" ]; \
  then git checkout develop && git pull && \
  git submodule update --init --recursive ; \
  else echo "on tag: $(git tag --contains)" ; \
fi
'
git diff --exit-code || \
(echo ERROR: Git changes detected. Submodules should either be tagged or on the latest version of develop. && exit 1)

# Make sure none of the submodules are ahead of upstream. Otherwise, the checks
# will fail when someone else pulls the repository
git submodule foreach \
    '
    if [ ${name} != "." ]; then
        git status | grep "Your branch is ahead" && \
        echo ERROR: Git detected ${name} is ahead of the remote. Please make sure all submodule changes have been pushed first. && exit 1
    fi
    echo ${name} is not ahead
    '
