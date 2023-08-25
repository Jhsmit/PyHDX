import time
from typing import Union, Optional
from pathlib import Path

import typer
from ipaddress import ip_address
from omegaconf import OmegaConf
from tqdm.auto import tqdm

from pyhdx.config import cfg
from pyhdx.datasets import DataVault

app = typer.Typer()


@app.command()
def serve(
    scheduler_address: Optional[str] = typer.Option(None, help="Address for dask scheduler to use"),
    config: Optional[Path] = typer.Option(
        None, exists=True, help="Optional PyHDX .yaml config file to use"
    ),
):
    """Launch the PyHDX web application"""

    from pyhdx.config import cfg
    from pyhdx.local_cluster import verify_cluster, default_cluster

    if config is not None:
        conf = OmegaConf.create(config.read_text())
        cfg.set_config(conf)

    if scheduler_address is not None:
        ip, port = scheduler_address.split(":")
        if not ip_address(ip):
            print("Invalid IP Address")
            return
        elif not 0 <= int(port) < 2**16:
            print("Invalid port, must be 0-65535")
            return
        cfg.cluster.scheduler_address = scheduler_address

    scheduler_address = cfg.cluster.scheduler_address
    if not verify_cluster(scheduler_address):
        # Start a new local cluster if none is found
        client = default_cluster()
        _, ip, port = client.scheduler_address.split(":")
        ip = ip.strip("/")
        scheduler_address = f"{ip}:{port}"
        print(f"Started new Dask LocalCluster at {scheduler_address}")

    # Start the PyHDX web application
    from pyhdx.web import serve as serve_pyhdx

    serve_pyhdx.run_apps()

    loop = True
    while loop:
        try:
            time.sleep(1)
        except KeyboardInterrupt:
            print("Interrupted")
            loop = False


datasets_app = typer.Typer(help="Manage HDX datasets")


@datasets_app.command()
def fetch(num: int = typer.Option(10, min=1, help="Maximum number of datasets to download")):
    """Update the datasets from the PyHDX repository"""
    vault = DataVault(cache_dir=cfg.database_dir)
    missing_datasets = set(vault.remote_index) - set(vault.datasets)

    if missing_datasets:
        todo = list(missing_datasets)[:num]
        for data_id in tqdm(todo):
            vault.fetch_dataset(data_id)


@datasets_app.command()
def clear():
    """Clear the local dataset cache"""
    vault = DataVault(cache_dir=cfg.database_dir)
    vault.clear_cache()


app.add_typer(datasets_app, name="datasets")


@app.callback()
def callback():
    pass


if __name__ == "__main__":
    app()
