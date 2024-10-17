####################################
Remote Development on MPCDF Clusters
####################################

Useful suggestions and instructions in the `MPCDF Technical Documentation <https://docs.mpcdf.mpg.de/>`_.




************************
MPCDF SSH Configurations
************************

In order to connect to the MPCDF clusters, you need to configure your SSH settings, and some SSH configuration tricks like setting up jump board make our life easier.
You can find the tricks from MPCDF document
`About gateway machine and tunneling <https://docs.mpcdf.mpg.de/faq/tricks.html#how-can-i-avoid-having-to-type-my-password-repeatedly-how-can-i-tunnel-through-the-gateway-machines>`_.

Here's an example SSH configuration for :code:`ADA` with :code:`mpcdf_gate2` as jump board without canonical matching:

.. code-block:: shell

    Host mpcdf_gate2
        HostName gate2.mpcdf.mpg.de
        User <your_username>
        ServerAliveInterval 30
        GSSAPIAuthentication yes
        GSSAPIDelegateCredentials yes
        ControlMaster auto
        ControlPersist 12h
        ControlPath ~/.ssh/master-%C

    Host a02
        HostName ada02.bc.rzg.mpg.de
        User <your_username>
        ProxyJump mpcdf_gate2
        ServerAliveInterval 30
        GSSAPIAuthentication yes
        GSSAPIDelegateCredentials yes
        ControlMaster auto
        ControlPersist 12h
        ControlPath ~/.ssh/master-%C

.. note::
   Replace `<your_username>` with your actual MPCDF username.




**************************************
VS Code + Remote Development Extension
**************************************

:code:`VS Code` is undoubtedly one of the most popular code editors among developers, featuring a wide range of extensions and plugins to enhance your development experience.
It is 100% free and semi-opensource, and available for all major operating systems.
You can install :code:`VS Code` from the official website `here <https://code.visualstudio.com/>`_.

.. Some introduction on how to install the VS Code Remote Development Extension.

One of the most useful extension packs for :code:`VS Code` is the `VS Code Remote Development <https://code.visualstudio.com/docs/remote/remote-overview>`_ extension suite.
It allows you to develop on a remote machine, container, or WSL, and you can use the same set of tools and extensions as you would on your local machine.




VS Code Setup on Remote Machine
###############################

Once you have installed :code:`VS Code`, you can install the `Remote Development extension <https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack>`_ in the `Extensions view <https://code.visualstudio.com/api/ux-guidelines/views>`_ by searching for `Remote Development` and install the packages by clicking on the `Install` button at the botton right corner of the plaquette.

Then the door to easy remote development is open for you: click on the monitor-like button, switch to :code:`Remotes (Tunnle/SSH)`, and then you should see all the SSH :code:`HostName` you have in your SSH config file.
But wait a second before you click the :code:`Connect in New Windows` - you don't want to type password and OTP over and over again, right?
I would suggest one always login to the target remote machine via a local terminal first, and the :code:`ControlMaster` settings helps you to keep the connection alive.
Then click the connection button, and we move on to the next section.

.. image:: https://code.visualstudio.com/assets/docs/remote/ssh/ssh-explorer-add-new.png
    :alt: VS Code Remote SSH Explorer
    :width: 600px
    :align: center

.. note::
    However, :code:`ControlMaster` is only available in Linux, so Windows users help yourself to some automation tools.




Suggested Extensions for Development
####################################

For Python development on MPCDF, we recommend the following VS Code extensions:

- `Python Extension Package`
    - The official Python extension for Visual Studio Code. It provides rich support for Python, including linting, IntelliSense, formatting, refactoring, debugging, unit tests, and Jupyter.
    - The `Python Environment Manager` is included in this extension suite - you can easily switch between different Python environments and manage your packages within the workspace.
- `Jupyter`
    - A method to proxy your SSH connection through an intermediate jump host.

Also we have some optional extensions for you to consider

- `Pylance`
    - It works alongside Python in Visual Studio Code to provide performant language support.
- `Ruff`
    - Even better Python linting and code formatting, but not yet popular.
- `Rainbow CSV`
    - Highlight CSV and TSV files in different colors, making it easier to read.
- `Trailing Spaces`
    - Develop your coding obsession by highlighting trailing spaces at the end of lines in bloody red.
- `Resource Monitor`
    - Monitor your system resources directly in VS Code, so you can kill your memory-hogging processes in time and IT team won't knock on your door.
- `Error Lens`
    - Directly show the error message inline, i.e. where the error occurs in your code, so you don't have to scroll up and down to find the error message.
- `Github Copilot` and `Github Copilot Chat`
    - The AI pair programming tool from Github, which can help you write code faster and more efficiently.
    - Most useful for Python, like doing repeated tasks and documentation writing (Yes! Copilot!).
- `Docs View`
    - Displays hover documentation in the sidebar or panel.




**************
Best Practices
**************

Rapid Development with Jupyter Notebooks
########################################

You can directly run Jupyter notebooks on the remote machine.
Just create a :code:`.ipynb` file and open it, then work as usual as in the jupyter notebooks, with more programming supports from the VS Code extensions.

Debugging
*********

First switch to a dark theme in VS Code to avoid attracting more bugs, then you can use the built-in debugger.
You can set breakpoints, step through your code, and inspect variables as debugging a usual python script.
The shortcut keybindings can be found in the :code:`Keyboard Shortcuts` settings by searching for `jupyter debug`.

Autoreload
**********

You must have seen the :code:`autoreload` magic command in Jupyter notebooks - it reloads the modules automatically before executing the code, so it is useful when you are developing a module and want to see the changes immediately.
You can find more `here <https://ipython.org/ipython-doc/3/config/extensions/autoreload.html>`_.

For example I have my own python package :code:`mypytools` and I am working on a file :code:`mypytools/utils.py`.
Then I can use the following commands in the Jupyter notebook to reload this python file automatically:

.. code-block:: python

    %load_ext autoreload
    %autoreload 1
    %aimport mypytools.utils
    from mypytools.utils import my_tool_func

and :code:`%autoreload 1` means "Reload all modules imported with %aimport every time before executing the Python code typed".




Other Tricks
************

You can toggle the line numbering by clicking on the blank area (switch to non-inputFocus status) and then do keyboard shortcut :code:`Shift+L`.

You can avoid super-long cell output by enabling `notebook.output.textLineLimit` in the settings, then you can have the outputs in boxes with scrollbars.


.. code-block:: python

   %reload_ext autoreload
   %autoreload 2

This will allow you to reload any updated modules during development.




Use GPU CUDA in Jupyter Notebooks
#################################

.. note::

    This is auto-generated content by GPT!

If you need to use GPU nodes on MPCDF, here’s how to request a GPU node with `salloc`:

.. code-block:: shell

   salloc --nodes=1 --gres=gpu:1 --time=01:00:00

Once the allocation is successful, start a Jupyter notebook on the allocated node:

.. code-block:: shell

   jupyter notebook --no-browser --port=8888

Then, on your local machine, you can connect using SSH tunneling:

.. code-block:: shell

   ssh -L 8888:localhost:8888 your_username@cluster.mpcdf.mpg.de

.. note::
   This allows you to utilize CUDA for GPU-accelerated computations in your Jupyter notebooks.


