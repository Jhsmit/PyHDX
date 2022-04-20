import time
from typing import Union, Optional
from pathlib import Path

import typer
from ipaddress import ip_address
import yaml


app = typer.Typer()


@app.command()
def serve(
    scheduler_address: Optional[str] = typer.Option(
        None, help="Address for dask scheduler to use"
    )
):
    """Launch the PyHDX web application"""

    from pyhdx.config import cfg
    from pyhdx.local_cluster import verify_cluster, default_cluster

    if scheduler_address is not None:
        ip, port = scheduler_address.split(":")
        if not ip_address(ip):
            print("Invalid IP Address")
            return
        elif not 0 <= int(port) < 2 ** 16:
            print("Invalid port, must be 0-65535")
            return
        cfg.set("cluster", "scheduler_address", scheduler_address)

    scheduler_address = cfg.get("cluster", "scheduler_address")
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


@app.command()
def process(
    jobfile: Path = typer.Argument(..., help="Path to .yaml jobfile"),
    cwd: Optional[Path] = typer.Option(None, help="Optional path to working directory"),
):
    """
    Process a HDX dataset according to a jobfile
    """

    from pyhdx.batch_processing import JobParser

    job_spec = yaml.safe_load(jobfile.read_text())
    parser = JobParser(job_spec, cwd=cwd)

    parser.execute()


if __name__ == "__main__":
    app()
