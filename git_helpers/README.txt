Post-commit and pre-commit are git hooks that help you increment the version number each time you
  commit. If you will be contributing to the GitHub repository, please copy these
  files to IGS_dada2_pipeline/.git/hooks and ensure they are executable.
OSX and Unix:
  cd PATH_TO_REMOTE/IGS_dada2_pipeline
  cp git_helpers/{post-commit,pre-commit} .git/hooks
  chmod 755 .git/hooks/{post-commit,pre-commit}

To make contributions that require review from other contributors, create a
  new branch on GitHub, clone that branch, and make changes on that branch.
  By initiating a pull request on GitHub, you will open up that branch for 
  review. When everyone is ready to switch to the new version of the pipeline,
  use GitHub to initiate a merge commit.

  * If you accidentally make changes in your clone of the wrong branch, you can 
  git stash save "name your changes in this message"
  git checkout new-branch
  git stash pop
