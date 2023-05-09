import time
from typing import Union, Optional
from pathlib import Path

import typer
from ipaddress import ip_address
import yaml
from omegaconf import OmegaConf

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


@app.callback()
def callback():
    pass


if __name__ == "__main__":
    app()
