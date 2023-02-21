"""Generate a pdf output with all peptide fits. Requires pdflatex"""
from pyhdx.output import FitReport
from pyhdx.fileIO import load_fitresult
from pathlib import Path
from concurrent import futures
import proplot as pplt


current_dir = Path().cwd()
fit_result = load_fitresult(current_dir / "output" / "SecB_tetramer_dimer_batch")

tmp_dir = Path(__file__).parent / "temp"
tmp_dir.mkdir(exist_ok=True)

if __name__ == "__main__":
    report = FitReport(fit_result, temp_dir=tmp_dir)
    report.add_standard_figure("peptide_coverage_figure")
    report.add_standard_figure("residue_time_scatter_figure")
    report.add_standard_figure("residue_scatter_figure")
    report.add_standard_figure("dG_scatter_figure", ncols=1, aspect=3)
    report.add_standard_figure("ddG_scatter_figure", ncols=1, reference=0)
    report.add_standard_figure(
        "linear_bars", cmap="viridis", norm=pplt.Norm("linear", 15e3, 35e3)
    )  # todo name from kwargs
    report.add_standard_figure("rainbowclouds")

    report.add_peptide_uptake_curves()

    executor = futures.ProcessPoolExecutor(max_workers=10)
    report.generate_figures(executor=executor)

    report.generate_latex()
    report.generate_pdf(current_dir / "output" / "fit_report")
