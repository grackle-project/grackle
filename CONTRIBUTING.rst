.. _contributing-code:

How to Develop Grackle
======================

Grackle is a community project!

We are very happy to accept patches, features, and bugfixes from any member of
the community!  Grackle is developed using mercurial, primarily because it
enables very easy and straightforward submission of changesets.  We're eager to
hear from you.

.. note:: If you already know how to use the `mercurial version control system
   <http://mercurial-scm.org>`_ and are comfortable with handling it yourself,
   the quickest way to contribute to Grackle is to `fork us on BitBucket
   <http://bitbucket.org/grackle/grackle/fork>`_, make your changes, push the
   changes to your fork and issue a `pull request
   <http://bitbucket.org/grackle/grackle/pull-requests>`_.  The rest of this
   document is just an explanation of how to do that.

Keep in touch, and happy hacking!

.. _open-issues:

Open Issues
-----------

If you're interested in participating in Grackle development, take a look at the
`issue tracker on bitbucket <https://bitbucket.org/grackle/issues>`_. If you are
encountering a bug that is not already tracked there, please `open a new issue
<https://bitbucket.org/grackle/grackle/issues/new>`_.

Submitting Changes
------------------

We provide a brief introduction to submitting changes here.  We encourage
contributions from any user.  While we do not discuss version control, mercurial
or the advanced usage of BitBucket in detail here, we do provide an outline of
how to submit changes and we are happy to provide further assistance or
guidance on the mailing list.

Licensing
+++++++++

Grackle is under the Enzo public license, a BSD-like license.

All contributed code must be BSD-compatible.  If you'd rather not license in
this manner, but still want to contribute, please consider creating an external
package, which we'll happily link to in the Grackle documentation.

How To Get The Source Code For Editing
++++++++++++++++++++++++++++++++++++++

Grackle is hosted on BitBucket, and you can see all of the Grackle repositories at
http://bitbucket.org/grackle/. In order to modify the source code for Grackle,
we ask that you make a "fork" of the main Grackle repository on bitbucket.  A
fork is simply an exact copy of the main repository (along with its history)
that you will now own and can make modifications as you please.  You can create
a personal fork by visiting the Grackle bitbucket webpage at
https://bitbucket.org/grackle/grackle/ .  After logging in, you should see an
option near the top right labeled "fork".  Click this option, and then click the
fork repository button on the subsequent page.  You now have a forked copy of
the Grackle repository for your own personal modification.

This forked copy exists on the bitbucket repository, so in order to access
it locally, follow the instructions at the top of that webpage for that
forked repository, namely run at a local command line:

.. code-block:: bash

   $ hg clone http://bitbucket.org/<USER>/<REPOSITORY_NAME>

This downloads that new forked repository to your local machine, so that you can
access it, read it, make modifications, etc.  It will put the repository in a
local directory of the same name as the repository in the current working
directory. You should also run the following command, to make sure you have the
most up-to-date version of Grackle checked out in your working directory.

.. code-block:: bash

   $ hg update default

You can see any past state of the code by using the hg log command.
For example, the following command would show you the last 5 changesets
(modifications to the code) that were submitted to that repository.

.. code-block:: bash

   $ cd <REPOSITORY_NAME>
   $ hg log -l 5

Using the revision specifier (the number or hash identifier next to each
changeset), you can update the local repository to any past state of the
code (a previous changeset or version) by executing the command:

.. code-block:: bash

   $ hg update revision_specifier

.. _mercurial-with-grackle:

How to Use Mercurial with Grackle
---------------------------------

If you're new to Mercurial, these three resources are pretty great for learning
the ins and outs:

* http://hginit.com/
* http://book.mercurial-scm.org
* http://mercurial-scm.org/
* http://mercurial-scm.org/wiki

The commands that are essential for using mercurial include:

* ``hg help`` which provides help for any mercurial command. For example, you
  can learn more about the ``log`` command by doing ``hg help log``. Other useful
  topics to use with ``hg help`` are ``hg help glossary``, ``hg help config``,
  ``hg help extensions``, and ``hg help revsets``.
* ``hg commit`` which commits changes in the working directory to the
  repository, creating a new "changeset object."
* ``hg add`` which adds a new file to be tracked by mercurial.  This does
  not change the working directory.
* ``hg pull`` which pulls (from an optional path specifier) changeset
  objects from a remote source.  The working directory is not modified.
* ``hg push`` which sends (to an optional path specifier) changeset objects
  to a remote source.  The working directory is not modified.
* ``hg log`` which shows a log of all changeset objects in the current
  repository.  Use ``-G`` to show a graph of changeset objects and their
  relationship.
* ``hg update`` which (with an optional "revision" specifier) updates the
  state of the working directory to match a changeset object in the
  repository.
* ``hg merge`` which combines two changesets to make a union of their lines
  of development.  This updates the working directory.

We are happy to asnswers questions about mercurial use on on the mailing list to
walk you through any troubles you might have.  Here are some general suggestions
for using mercurial:

* Named branches are to be avoided.  Try using bookmarks (``see hg help
  bookmark``) to track work.  (`More info about bookmarks is available on the
  mercurial wiki <http://mercurial-scm.org/wiki/Bookmarks>`_)
* Make sure you set a username in your ``~/.hgrc`` before you commit any
  changes!  All of the tutorials above will describe how to do this as one of
  the very first steps.
* Please avoid deleting your Grackle forks, as that deletes the pull request
  discussion from process from BitBucket's website, even if your pull request
  is merged.
* You should only need one fork. See :ref:`sharing-changes` for a description of
  the basic workflow.

.. _sharing-changes:

Making and Sharing Changes
--------------------------

The simplest way to submit changes to Grackle is to do the following:

* Build Grackle from the mercurial repository
* Navigate to the root of the Grackle repository
* Make some changes and commit them
* Fork the `Grackle repository on BitBucket
  <https://bitbucket.org/grackle/grackle>`_
* Push the changesets to your fork
* Issue a pull request.

Here's a more detailed flowchart of how to submit changes.

#. Edit the source file you are interested in and test your changes.
#. Fork Grackle on BitBucket.  (This step only has to be done once.)  You can do
   this at: https://bitbucket.org/grackle/grackle/fork.  Call this repository
   grackle.
#. Create a bookmark to track your work. For example: ``hg bookmark
   my-first-pull-request``
#. Commit these changes, using ``hg commit``.  This can take an argument
   which is a series of filenames, if you have some changes you do not want
   to commit.
#. Remember that this is a large development effort and to keep the code
   accessible to everyone, good documentation is a must.  Add in source code
   comments for what you are doing.  Add in docstrings
   if you are adding a new function or class or keyword to a function.
   Add documentation to the appropriate section of the online docs so that
   people other than yourself know how to use your new code.
#. If your changes include new functionality or cover an untested area of the
   code, add a test. Commit these changes as well.
#. Push your changes to your new fork using the command::

      hg push -B my-first-pull-request https://bitbucket.org/YourUsername/grackle/

   Where you should substitute the name of the bookmark you are working on for
   ``my-first-pull-request``. If you end up doing considerable development, you
   can set an alias in the file ``.hg/hgrc`` to point to this path.

   .. note::
     Note that the above approach uses HTTPS as the transfer protocol
     between your machine and BitBucket.  If you prefer to use SSH - or
     perhaps you're behind a proxy that doesn't play well with SSL via
     HTTPS - you may want to set up an `SSH key`_ on BitBucket.  Then, you use
     the syntax ``ssh://hg@bitbucket.org/YourUsername/grackle``, or equivalent,
     in place of ``https://bitbucket.org/YourUsername/grackle`` in Mercurial
     commands. For consistency, all commands we list in this document use the
     HTTPS protocol.

     .. _SSH key: https://confluence.atlassian.com/display/BITBUCKET/Set+up+SSH+for+Mercurial

#. Issue a pull request at
   https://bitbucket.org/YourUsername/grackle/pull-request/new
   A pull request is an automated way of asking people to review and accept the
   modifications you have made to your personal version of the code.

During the course of your pull request you may be asked to make changes.  These
changes may be related to style issues, correctness issues, or requesting
tests.  The process for responding to pull request code review is relatively
straightforward.

#. Make requested changes, or leave a comment on the pull request page on
   Bitbucket indicating why you don't think they should be made.
#. Commit those changes to your local repository.
#. Push the changes to your fork:

      hg push https://bitbucket.org/YourUsername/grackle/

#. Your pull request will be automatically updated.
