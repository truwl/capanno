First release checklist
- [ ] Make sure a `gh-pages` branch is setup.
  - [ ] Contains the `index.html`.
  - [ ] The values in the `_config.yml` are correct.
- [ ] Make sure Zenodo is enabled for the repo.
  - [ ] Switch on at https://zenodo.org/account/settings/github/.
  - [ ] Make sure the zenodo record ID is set in the `_config.yml` on the `gh-pages` branch. This ID can be found by clicking on the DOI badge displayed under the relevant repo at https://zenodo.org/account/settings/github/. The number in the urls that get displayed is the number that should be put into the config.
- [ ] Check that the repository tagline has the appropriate category at the end, eg. ` Category: Multi-Sample`.
- [ ] Follow the general checklist.
- [ ] Update the biowdl homepage, by pushing on empty commit to the master branch of biowdl/biowdl.github.io.

General release checklist
- [ ] Check outstanding issues on JIRA and Github
- [ ] Update all submodules to latest master with: `git submodule foreach "git checkout master;git pull; git submodule foreach --recursive 'git fetch'; git submodule update --init --recursive"`
- [ ] check all submodules are tagged correctly with `git submodule`
- [ ] run tests to confirm to be released version works.
- [ ] Generate inputs overview using wdl-aid:
  `wdl-aid --strict -t scripts/docs_template.md.j2 pipeline.wdl > docs/inputs.md`
- [ ] Publish documentation (`updateDocs.sh`) from `develop` branch
  - [ ] Copy docs folder to `gh-pages` branch
  - [ ] Overwrite existing develop folder with docs folder on `gh-pages`
  - [ ] Push changes to `gh-pages branch`
- [ ] Check [latest documentation
](https://biowdl.github.io/) looks fine
- [ ] Change current development version in `CHANGELOG.md` to stable version.
- [ ] Run the release script `release.sh`
  - [ ] Check all submodules are tagged
  - [ ] Merge the develop branch into `master`
  - [ ] Created an annotated tag with the stable version number. Include changes 
    from changelog.md.
  - [ ] Confirm or set stable version to be used for tagging
  - [ ] Push tag to remote.
  - [ ] Merge `master` branch back into `develop`.
  - [ ] Add updated version number to develop
- [ ] Publish documentation (`updateDocs.sh`) from `master` branch
  - [ ] Copy docs folder to `gh-pages` branch
  - [ ] Rename docs to new stable version on `gh-pages`
  - [ ] Set latest version to new version
  - [ ] Push changes to `gh-pages branch`
- [ ] Create a new release from the pushed tag on github
- [ ] Prepare the repo for packaging by `git checkout master && git submodule update --init --recursive`
  - [ ] Package the wdl files with `wdl-packager --reproducible -a LICENSE 
      -a dockerImages.yml <WDL_FILE>`
  - [ ] Add the package(s) to the github release. Also add the original WDL file
      as `<pipeline>_<version>.wdl` following the same naming as the package.
      This alllows for usage of wdl and imports zip with cromwell without 
      requiring the user to extract the package.

  
  
