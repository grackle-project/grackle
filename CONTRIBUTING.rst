.. _contributing-code:

How to Develop Grackle
======================

Grackle is a community project!

We are very happy to accept patches, features, and bugfixes from any member of
the community!  Grackle is developed using Git, primarily because of how well
it enables open-source, community contribution. We're eager to hear from you.

.. note:: If you are already familiar with Git and `GitHub <https://github.com>`_,
   the best way to contribute is to fork the `main Grackle repository
   <https://github.com/grackle-project/grackle>`__, make your changes, push them
   to your fork, and issue a pull request. The rest of this document is just an
   explanation of how to do that.

Keep in touch, and happy hacking!

.. _open-issues:

Open Issues
-----------

If you're interested in participating in Grackle development, take a look at the
`issue tracker on GitHub <https://github.com/grackle-project/grackle/issues>`_.
If you are encountering a bug that is not already tracked there, please `open a
new issue <https://github.com/grackle-project/grackle/issues/new>`__.

Contributing to Grackle with Git and Github
-------------------------------------------

We provide a brief introduction to submitting changes here.  We encourage
contributions from any user. If you are new to Git and/or GitHub, there are
excellent guides available at `guides.github.com <https://guides.github.com/>`_,
specifically the `Git Handbook
<https://guides.github.com/introduction/git-handbook/>`__, and the `GitHub
Hello World <https://guides.github.com/activities/hello-world/>`__. We are also
happy to provide guidance on the mailing list or in our slack channel.

Licensing
+++++++++

Grackle is under the Enzo public license, a BSD-like license.

All contributed code must be BSD-compatible.  If you'd rather not license in
this manner, but still want to contribute, please consider creating an external
package, which we'll happily link to in the Grackle documentation.

How To Get The Source Code For Editing
++++++++++++++++++++++++++++++++++++++

Grackle is hosted on GitHub, and you can see all of the Grackle repositories at
https://github.com/grackle-project/. In order to modify the source code for Grackle,
we ask that you make a "fork" of the main Grackle repository on GitHub.  A
fork is simply an exact copy of the main repository (along with its history)
that you will now own and can make modifications as you please.  You can create
a personal fork by visiting the Grackle GitHub webpage at
https://github.com/grackle-project/grackle/.  After logging in, you should see an
option near the top right labeled "fork".  You now have a forked copy of
the Grackle repository for your own personal modification.

This forked copy exists on GitHub under your username, so in order to access
it locally, follow the instructions at the top of that webpage for that
forked repository:

.. code-block:: bash

   $ git clone http://bitbucket.org/<USER>/<REPOSITORY_NAME>

This downloads that new forked repository to your local machine, so that you can
access it, read it, make modifications, etc.  It will put the repository in a
local directory of the same name as the repository in the current working
directory.

.. code-block:: bash

   $ cd grackle

Verify that you are on the main branch of Grackle by running:

.. code-block:: bash

   $ git branch

If you're not on the main branch, you can get to it with:

.. code-block:: bash

   $ git checkout main

You can see any past state of the code by using the git log command.
For example, the following command would show you the last 5 revisions
(modifications to the code) that were submitted to that repository.

.. code-block:: bash

   $ git log -n 5

Using the revision specifier (the number or hash identifier next to each
changeset), you can update the local repository to any past state of the
code (a previous changeset or version) by executing the command:

.. code-block:: bash

   $ git checkout revision_specifier

.. _sharing-changes:

Making and Sharing Changes
--------------------------

The simplest way to submit changes to Grackle is to do the following:

#. Fork the main repository.
#. Clone your fork.
#. Make some changes and commit them.
#. Push the changesets to your fork.
#. Issue a pull request.

Here's a more detailed flowchart of how to submit changes.

#. Fork Grackle on GitHub.  (This step only has to be done once.)  You can do
   this by clicking on the **fork** button in the top-right corner of `the main
   repository <https://github.com/grackle-project/grackle>`__.
#. Create a new branch in which to make your changes by doing ``git
   checkout -b <new branch name>``. This will make it easy to move back and
   forth between the main branch of the code and your changes.
#. Edit the source file you are interested in and test your changes.
#. Use ``git add <files>`` to stage files to be committed.
#. Commit your changes with ``git commit``. This will open a text editor so you
   can write a commit message. To add your message inline, do
   ``git commit -m "<commit message>"``. You can list specific file to be
   committed.
#. Remember that this is a large development effort and to keep the code
   accessible to everyone, good documentation is a must.  Add in source code
   comments for what you are doing.  Add documentation to the appropriate
   section of the online docs so that people other than yourself know how
   to use your new code.
#. If your changes include new functionality or cover an untested area of the
   code, add a test. Commit these changes as well.
#. Push your changes to your new fork using the command::

      $ git push origin <branch name>

   .. note::
     Note that the above approach uses HTTPS as the transfer protocol
     between your machine and GitHub.  If you prefer to use SSH - or
     perhaps you're behind a proxy that doesn't play well with SSL via
     HTTPS - you may want to set up an `SSH key
     <https://help.github.com/articles/connecting-to-github-with-ssh/>`__
     on GitHub.  Then, you use
     the syntax ``ssh://git@github.com/<USER>/grackle``, or equivalent, in
     place of ``https://github.com/<USER>/grackle`` in git commands.
     For consistency, all commands we list in this document will use the HTTPS
     protocol.

#. Issue a pull request by going to the main repository and clicking on the
   green button that says **Compare & pull request**. This will open up a page
   that will allow you to enter a description of the changes to be merged. Once
   submitted, a series of automated tests will run and their status will be
   reported on the pull request page.

During the course of your pull request you may be asked to make changes.  These
changes may be related to style issues, correctness issues, or requesting
tests.  The process for responding to pull request code review is relatively
straightforward.

#. Make requested changes, or leave a comment on the pull request page on
   GitHub indicating why you don't think they should be made.
#. Commit those changes to your local repository.
#. Push the changes to your fork::

      $ git push origin <branch name>

#. Your pull request will be automatically updated.

Once your pull request has been accepted, you can safely delete your
branch::

      $ git branch --delete <branch name>

Updating Your Branch
++++++++++++++++++++

If your branch or pull request has been open for some time, it may be useful
to keep it up to date with the latest changes from the main repository. This
can be done by `rebasing your changes <https://git-scm.com/docs/git-rebase>`__.
Before doing this, you will need to be able to pull the latest changes from
the main repository.

#. Add the main repository as a remote::

      $ git remote add grackle https://github.com/grackle-project/grackle

   You can verify that it has been added by doing ``git remote -v``. This
   only needs to be done once.

#. Go back to the main branch and pull the changes::

      $ git checkout main
      $ git pull grackle main

#. Return to your branch and rebase your changes onto the head of the main
   branch::

      $ git checkout <branch name>
      $ git rebase main

This should go smoothly unless changes have been made to the same lines in
the source, in which case you will need to fix conflicts. After rebasing,
you will get an error when trying to push your branch to your fork. This is
because you have changed the order of commits and git does not like that.
In this case, you will need to add "-f" to your push command to force
the changes to be accepted.::

      $ git push -f origin <branch name>

Have fun!