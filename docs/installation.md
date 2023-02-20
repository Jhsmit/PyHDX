# Installation

## Stable release

Installation of the latest stable release version.

With `conda`:

```bash
conda install -c conda-forge pyhdx
```

With `pip`:

```bash
pip install pyhdx
```

To install with web application:

```bash
pip install pyhdx[web]
```

To install with pdf output:
```bash
pip install pyhdx[pdf]
```

Installation via `conda` automatically installs all extras. 

## Running the web server


PyHDX web application can be launched from the command line using  the `pyhdx` cli command with below options,

To run PyHDX server using default settings on your local server:

```bash
pyhdx serve
```

To run PyHDX server using the IP address and port number of your dask cluster:

```bash
pyhdx serve --scheduler_address <ip>:<port>
```


If no dask cluster is found at the specified address, a LocalCluster will be started (on localhost) using the
specified port number.

To start a dask cluster separately, open another terminal tab and run:

```bash
python local_cluster.py
```

This will start a Dask cluster on the scheduler address as specified in the PyHDX config.


## Install from source


Create a new conda environment:

```bash
conda create --name py39_pyhdx python=3.9
conda activate py39_pyhdx
```

Clone the GitHub repository:
```bash
git clone https://github.com/Jhsmit/PyHDX
cd PyHDX
```

Dependencies can then be installed with `poetry`

```bash
    $ poetry install --all-extras
```

Use `--all-extras` if you plan to use the web interface. 


### Running from source

To run the web application:

```bash
python pyhdx/web/serve.py
```

This runs the pyhx web application without a Dask cluster to submit jobs to, so
submitting a fitting job will give an error.

To start a dask cluster separately, open another terminal tab and run:

```bash
python pyhdx/local_cluster.py
```

## Configuration

A configuration file is located in the `.pyhdx` folder in the user home directory. This file
is used by default and can be edited to change PyHDX default settings.

Alternatively, users can create additional `.yaml` configuration files in this directory, after
which the scripts ``local_cluster.py`` and ``serve.py`` prompt the user for which file to use.

The section ``server`` configures the panel server settings. In this section the additional keys
``port`` and ``websocket_origin`` can be added, which are passed to ``panel.serve``. See the panel 
[Deploy and Export](https://panel.holoviz.org/user_guide/Deploy_and_Export.html) deploy section for 
more information.

